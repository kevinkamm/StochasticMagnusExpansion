clear all;gpuDevice(1);close all
mExpDevice='cpuReshape';
%% Parameters
% Gloabal
param.t0=0;
param.T=1;
param.d=2;
[A,B]=AB_const_coeff_initialize(param,'Example','default');
moments=[1,2,3];
% Z=2.576;% 99% confidence interval
Z=1.96; % 95% confidence interval 
% Reference method
deltaEulerRef=10^(-4);
sampleEulerRef = 10^4;

% Convergence test
deltaStart = 10^(-3);
deltaEnd = 10^(-4);
deltaSteps = 2;
deltaConv = linspace(deltaStart,deltaEnd,deltaSteps);
sampleStart = 10^3;
sampleEnd = 10^4;
sampleSteps = 2;
sampleConv = [10^3,10^4];%linspace(sampleStart,sampleEnd,sampleSteps);
% Result variables
ctimeEulerRef=0;
ctimeEulerRefTotal=zeros(length(moments),1);
momentEulerRef=zeros(param.d,param.d,length(moments));
ctimeEulerConv=zeros(sampleSteps,deltaSteps);
ctimeEulerConvTotal=zeros(length(moments),sampleSteps,deltaSteps);
momentEulerConv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
rmseEulerConv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
ctimeM1ConvLog=zeros(sampleSteps,deltaSteps);
ctimeM1ConvExp=zeros(sampleSteps,deltaSteps);
ctimeM1ConvTotal=zeros(length(moments),sampleSteps,deltaSteps);
momentM1Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
rmseM1Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
ctimeM2ConvLog=zeros(sampleSteps,deltaSteps);
ctimeM2ConvExp=zeros(sampleSteps,deltaSteps);
ctimeM2ConvTotal=zeros(length(moments),sampleSteps,deltaSteps);
momentM2Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
rmseM2Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
ctimeM3ConvLog=zeros(sampleSteps,deltaSteps);
ctimeM3ConvExp=zeros(sampleSteps,deltaSteps);
ctimeM3ConvTotal=zeros(length(moments),sampleSteps,deltaSteps);
momentM3Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
rmseM3Conv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
% Save paths
root=[pwd,'\','Moment','\','Matlab'];

%% Computing reference method
disp('Computing reference method')
% initialize
param.N=ceil(deltaEulerRef^(-1)*param.T);
param.M=sampleEulerRef;
[time,BM]=AB_const_param_initialize(param);

% Calculate euler
disp('Calculate euler')
ticEulerRef=tic;
XeulerRef=AB_const_euler(time.dt,...
                    BM.dW,...
                    A,...
                    B,...
                    'Comp Device','gpu', 'disp', 0);
ctimeEulerRef=toc(ticEulerRef);               
fprintf('Elapsed time %g\n',ctimeEulerRef);
% Calculate moments
disp('Calculate moments')
for mm=1:1:length(moments)
    morder=moments(mm);
    ticOrder=tic;
    momentEulerRef(:,:,mm)=mean(XeulerRef(:,:,end,:).^morder,4);
    ctimeEulerRefTotal(mm)=ctimeEulerRef+toc(ticOrder);
end
clear XeulerRef
%% Save moments and ctimes
folderName=sprintf('momentRef_T%1.3f_d%d_N%d_M%d',...
            param.T,param.d,param.N,param.M);
folderPath=[root,'\',folderName];
if exist(folderPath)~=7
    mkdir(folderPath);
end
file=sprintf('momentEuler');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'momentEulerRef')
file=sprintf('ctimeEuler');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeEulerRef')
file=sprintf('ctimeEulerTotal');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeEulerRefTotal')

%% Computing convergence test
for isample=1:1:length(sampleConv)
    for idelta=1:1:length(deltaConv)
        gpuDevice(1);
        fprintf('Runs convergence: %d/%d samples, %d/%d deltas\n',...
            isample,length(sampleConv),idelta,length(deltaConv))
        % initialize
        param.N=ceil(deltaConv(idelta)^(-1)*param.T);
        param.M=sampleConv(isample);
        [time,BM]=AB_const_param_initialize(param);
        dt=time.dt;
        dW=BM.dW;
        % Calculate euler
        disp('Calculate euler')
        ticEuler=tic;
        X=AB_const_euler(dt,...
                            dW,...
                            A,...
                            B,...
                            'Comp Device','gpu');
        ctimeEulerConv(isample,idelta)=toc(ticEuler);
        fprintf('Elapsed time %g\n',ctimeEulerConv(isample,idelta));
        % Calculate moments for euler
        disp('Calculate moments for euler')
        for mm=1:1:length(moments)
            morder=moments(mm);
            ticOrder=tic;
            momentEulerConv(:,:,mm,isample,idelta)=mean(X(:,:,end,:).^morder,4);
            ctimeEulerConvTotal(mm,isample,idelta)=ctimeEulerConv(isample,idelta)+toc(ticOrder);
            % Calculate rmse for euler
            rmseEulerConv(:,:,mm,isample,idelta)=...
                sqrt(var(X(:,:,end,:).^morder-momentEulerRef(:,:,mm),0,4)./size(X,4));
        end
%         gpuDevice(1);
        N=ceil((10^(-2))^(-1)*param.T);
        ind_t=[1:1:N];
        ind_t(2:1:end)=ind_t(1:1:end-1).*floor((param.N-1)/(N-1))+1;

        t=linspace(param.t0,param.T,N);
        fprintf('diff of t %g\n',sum(abs(reshape(time.t(1,1,ind_t,1),[],1)-reshape(t,[],1)),'all'))
        dt=(param.T-param.t0)/(N-1);
        W=BM.W(:,:,ind_t,:);
        % Calculate Magnus order 1
        disp('Calculate Magnus order 1')
        % Calculate log of Magnus order 1
        disp('Calculating log\n');
        ticMagnus=tic;
        O=AB_const_magnus(...
                time.t(1,1,ind_t,1),...
                W,...
                A,...
                B,...
                1,...
                'Comp Device','gpu');
        ctimeM1ConvLog(isample,idelta)=toc(ticMagnus);
        fprintf('Elapsed time %g\n',ctimeM1ConvLog(isample,idelta));
        % Calculate exp of Magnus order 1
        disp('Calculating exp\n');
        ticExp=tic;
        X=m_exp(O(:,:,end,:),...
            'type','matlab',...
            'Comp Device',mExpDevice);
        ctimeM1ConvExp(isample,idelta)=toc(ticExp);
        fprintf('Elapsed time %g\n',ctimeM1ConvExp(isample,idelta));
        % Calculate moments for Magnus order 1
        disp('Calculate moments for Magnus order 1')
        for mm=1:1:length(moments)
            morder=moments(mm);
            ticOrder=tic;
            momentM1Conv(:,:,mm,isample,idelta)=mean(X(:,:,end,:).^morder,4);
            ctimeM1ConvTotal(mm,isample,idelta)=ctimeM1ConvLog(isample,idelta)+...
                ctimeM1ConvExp(isample,idelta)+toc(ticOrder);
            % Calculate rmse for M1 wrt reference method
            rmseM1Conv(:,:,mm,isample,idelta)=...
                sqrt(var(X(:,:,end,:).^morder-momentEulerRef(:,:,mm),0,4)./size(X,4));
        end
%         tic;
%         toc;
%         gpuDevice(1);
        % Calculate Magnus order 2
        disp('Calculate Magnus order 2')
        % Calculate log of Magnus order 2
        disp('Calculating log\n');
        ticMagnus=tic;
        O=AB_const_magnus(...
                time.t(1,1,ind_t,1),...
                W,...
                A,...
                B,...
                2,...
                'Comp Device','gpu');
        ctimeM2ConvLog(isample,idelta)=toc(ticMagnus);
        fprintf('Elapsed time %g\n',ctimeM2ConvLog(isample,idelta));
        % Calculate exp of Magnus order 1
        disp('Calculating exp\n');
        ticExp=tic;
        X=m_exp(O(:,:,end,:),...
            'type','matlab',...
            'Comp Device',mExpDevice);
        ctimeM2ConvExp(isample,idelta)=toc(ticExp);
        fprintf('Elapsed time %g\n',ctimeM2ConvExp(isample,idelta));
%         tic;
%         toc;
        % Calculate moments for Magnus order 2
        disp('Calculate moments for Magnus order 2')
        for mm=1:1:length(moments)
            morder=moments(mm);
            ticOrder=tic;
            momentM2Conv(:,:,mm,isample,idelta)=mean(X(:,:,end,:).^morder,4);
            ctimeM2ConvTotal(mm,isample,idelta)=ctimeM2ConvLog(isample,idelta)+...
                ctimeM2ConvExp(isample,idelta)+toc(ticOrder);
            % Calculate rmse for M2 wrt reference method
            rmseM2Conv(:,:,mm,isample,idelta)=...
                sqrt(var(X(:,:,end,:).^morder-momentEulerRef(:,:,mm),0,4)./size(X,4));
        end
%         gpuDevice(1);
%         tic;
%         toc;
        % Calculate Magnus order 3
        disp('Calculate Magnus order 3')
        % Calculate log of Magnus order 3
        disp('Calculating log\n');
        ticMagnus=tic;
        O=AB_const_magnus(...
                time.t(1,1,ind_t,1),...
                W,...
                A,...
                B,...
                3,...
                'Comp Device','gpu');
        ctimeM3ConvLog(isample,idelta)=toc(ticMagnus);
        fprintf('Elapsed time %g\n',ctimeM3ConvLog(isample,idelta));
        % Calculate exp of Magnus order 3
        disp('Calculating exp\n');
        ticExp=tic;
        X=m_exp(O(:,:,end,:),...
            'type','matlab',...
            'Comp Device',mExpDevice);
        ctimeM3ConvExp(isample,idelta)=toc(ticExp);
        fprintf('Elapsed time %g\n',ctimeM3ConvExp(isample,idelta));
        % Calculate moments for Magnus order 3
        disp('Calculate moments for Magnus order 3')
        for mm=1:1:length(moments)
            morder=moments(mm);
            ticOrder=tic;
            momentM3Conv(:,:,mm,isample,idelta)=mean(X(:,:,end,:).^morder,4);
            ctimeM3ConvTotal(mm,isample,idelta)=ctimeM3ConvLog(isample,idelta)+...
                ctimeM3ConvExp(isample,idelta)+toc(ticOrder);
            % Calculate rmse for M3 wrt reference method
            rmseM3Conv(:,:,mm,isample,idelta)=...
                sqrt(...
                    var(X(:,:,end,:).^morder-momentEulerRef(:,:,mm),0,4)./...
                    size(X,4));
        end
    end
end
clear O X
%% Save moments and ctimes
% Euler
folderName=sprintf('momentConv_T%1.3f_d%d_N%d_M%d',...
            param.T,param.d,param.N,param.M);
folderPath=[root,'\',folderName];
if exist(folderPath)~=7
    mkdir(folderPath);
end
file=sprintf('momentEuler');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'momentEulerConv')
file=sprintf('rmseEuler');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'rmseEulerConv')
file=sprintf('ctimeEuler');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeEulerConv')
file=sprintf('ctimeEulerTotal');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeEulerConvTotal')
% Magnus order 1
% folderName=sprintf('momentConv_T%1.3f_d%d_N%d_M%d',...
%             param.T,param.d,param.N,param.M);
folderPath=[root,'\',folderName];
if exist(folderPath)~=7
    mkdir(folderPath);
end
file=sprintf('momentM1');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'momentM1Conv')
file=sprintf('rmseM1');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'rmseM1Conv')
file=sprintf('ctimeM1Log');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM1ConvLog')
file=sprintf('ctimeM1Exp');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM1ConvExp')
file=sprintf('ctimeM1Total');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM1ConvTotal')
% Magnus order 2
% folderName=sprintf('momentConv_T%1.3f_d%d_N%d_M%d',...
%             param.T,param.d,param.N,param.M);
folderPath=[root,'\',folderName];
if exist(folderPath)~=7
    mkdir(folderPath);
end
file=sprintf('momentM2');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'momentM2Conv')
file=sprintf('rmseM2');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'rmseM2Conv')
file=sprintf('ctimeM2Log');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM2ConvLog')
file=sprintf('ctimeM2Exp');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM2ConvExp')
file=sprintf('ctimeM2Total');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM2ConvTotal')
% Magnus order 3
% folderName=sprintf('momentConv_T%1.3f_d%d_N%d_M%d',...
%             param.T,param.d,param.N,param.M);
folderPath=[root,'\',folderName];
if exist(folderPath)~=7
    mkdir(folderPath);
end
file=sprintf('momentM3');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'momentM3Conv')
file=sprintf('rmseM3');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'rmseM3Conv')
file=sprintf('ctimeM3Log');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM3ConvLog')
file=sprintf('ctimeM3Exp');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM3ConvExp')
file=sprintf('ctimeM3Total');
savePath=[folderPath,'\',file];
if exist(savePath)
    delete(savePath);
end
save(savePath,'ctimeM3ConvTotal')
%% Plots
% rsme  
for mm=1:1:length(moments)
    morder=moments(mm);
    fig = figure('units','normalized',...
              'outerposition',[0 0 1 1]); hold on;
    fontsize=22;
    linewidth=1;
    set(gca,'FontSize',fontsize)
    set(fig,'defaultlinelinewidth',linewidth)
    set(fig,'defaultaxeslinewidth',linewidth)
    set(fig,'defaultpatchlinewidth',linewidth)
    set(fig,'defaultAxesFontSize',fontsize)
    for row=1:1:param.d
        for col=1:1:param.d
            nexttile;hold on;
            for isample=1:1:length(sampleConv)
                color.euler = [1 0 0] .* isample ./ length(sampleConv);
                color.m1 = [0 0 1] .* isample ./ length(sampleConv);
                color.m2 = [1 1 0] .* isample ./ length(sampleConv);
                color.m3 = [0 1 0] .* isample ./ length(sampleConv);
                plot(deltaConv,squeeze(rmseEulerConv(row,col,mm,isample,:)),...
                    'LineStyle','-','Color',color.euler);
%                 plot(deltaConv,squeeze(rmseM1Conv(row,col,mm,isample,:)),...
%                     'LineStyle','-','Color',color.m1);
                plot(deltaConv,squeeze(rmseM2Conv(row,col,mm,isample,:)),...
                    'LineStyle','-','Color',color.m2);
                plot(deltaConv,squeeze(rmseM3Conv(row,col,mm,isample,:)),...
                    'LineStyle','-','Color',color.m3);
            end
            xlabel('$\Delta$', 'fontweight', 'bold','Interpreter','Latex')
            ylabel(sprintf('RMSE for $X^{%d%d}$',row,col),...
                'fontweight', 'bold','Interpreter','Latex')
        end
    end
end
for mm=1:1:length(moments)
    morder=moments(mm);
    fig = figure('units','normalized',...
              'outerposition',[0 0 1 1]); hold on;
    fontsize=22;
    linewidth=1;
    set(gca,'FontSize',fontsize)
    set(fig,'defaultlinelinewidth',linewidth)
    set(fig,'defaultaxeslinewidth',linewidth)
    set(fig,'defaultpatchlinewidth',linewidth)
    set(fig,'defaultAxesFontSize',fontsize)
    for row=1:1:param.d
        for col=1:1:param.d
            nexttile;hold on;
            [mDelta,mSamples]=meshgrid(deltaConv,sampleConv);
                color.euler = [1 0 0];
                color.m1 = [0 0 1];
                color.m2 = [1 1 0];
                color.m3 = [0 1 0];
                surf(mDelta,mSamples,squeeze(rmseEulerConv(row,col,mm,:,:)),...
                    'LineStyle','none','FaceColor',color.euler,'FaceAlpha',.1);
%                 surf(mDelta,mSamples,squeeze(rmseM1Conv(row,col,mm,:,:)),...
%                     'LineStyle','none','FaceColor',color.m1,'FaceAlpha',.1);
                surf(mDelta,mSamples,squeeze(rmseM2Conv(row,col,mm,:,:)),...
                    'LineStyle','none','FaceColor',color.m2,'FaceAlpha',.1);
                surf(mDelta,mSamples,squeeze(rmseM3Conv(row,col,mm,:,:)),...
                    'LineStyle','none','FaceColor',color.m3,'FaceAlpha',.1);
            view(3);
            xlabel('Time step size', 'fontweight', 'bold','Interpreter','Latex')
            ylabel('Samples', 'fontweight', 'bold','Interpreter','Latex')
            zlabel(sprintf('RMSE for $X^{%d%d}$',row,col),...
                'fontweight', 'bold','Interpreter','Latex')
        end
    end
end
% computational time
for mm=1:1:length(moments)
    morder=moments(mm);
    fig = figure('units','normalized',...
              'outerposition',[0 0 1 1]); hold on;
    fontsize=22;
    linewidth=1;
    set(gca,'FontSize',fontsize)
    set(fig,'defaultlinelinewidth',linewidth)
    set(fig,'defaultaxeslinewidth',linewidth)
    set(fig,'defaultpatchlinewidth',linewidth)
    set(fig,'defaultAxesFontSize',fontsize)
    for isample=1:1:length(sampleConv)
        color.euler = [1 0 0] .* isample ./ length(sampleConv);
        color.m1 = [0 0 1] .* isample ./ length(sampleConv);
        color.m2 = [1 1 0] .* isample ./ length(sampleConv);
        color.m3 = [0 1 0] .* isample ./ length(sampleConv);
        plot(deltaConv,squeeze(ctimeEulerConvTotal(mm,isample,:)),...
            'LineStyle','-','Color',color.euler);
        plot(deltaConv,squeeze(ctimeM1ConvTotal(mm,isample,:)),...
            'LineStyle','-','Color',color.m1);
        plot(deltaConv,squeeze(ctimeM2ConvTotal(mm,isample,:)),...
            'LineStyle','-','Color',color.m2);
        plot(deltaConv,squeeze(ctimeM3ConvTotal(mm,isample,:)),...
            'LineStyle','-','Color',color.m3);
    end
    xlabel('$\Delta$', 'fontweight', 'bold','Interpreter','Latex')
    ylabel(sprintf('Comp Times  for $X^{%d%d}$',row,col),...
        'fontweight', 'bold','Interpreter','Latex')
end
for mm=1:1:length(moments)
    morder=moments(mm);
    fig = figure('units','normalized',...
              'outerposition',[0 0 1 1]); hold on;
    fontsize=22;
    linewidth=1;
    set(gca,'FontSize',fontsize)
    set(fig,'defaultlinelinewidth',linewidth)
    set(fig,'defaultaxeslinewidth',linewidth)
    set(fig,'defaultpatchlinewidth',linewidth)
    set(fig,'defaultAxesFontSize',fontsize)
    [mDelta,mSamples]=meshgrid(deltaConv,sampleConv);
        color.euler = [1 0 0];
        color.m1 = [0 0 1];
        color.m2 = [1 1 0];
        color.m3 = [0 1 0];
        surf(mDelta,mSamples,squeeze(ctimeEulerConvTotal(mm,:,:)),...
            'LineStyle','none','FaceColor',color.euler,'FaceAlpha',.1);
        surf(mDelta,mSamples,squeeze(ctimeM1ConvTotal(mm,:,:)),...
            'LineStyle','none','FaceColor',color.m1,'FaceAlpha',.1);
        surf(mDelta,mSamples,squeeze(ctimeM2ConvTotal(mm,:,:)),...
            'LineStyle','none','FaceColor',color.m2,'FaceAlpha',.1);
        surf(mDelta,mSamples,squeeze(ctimeM3ConvTotal(mm,:,:)),...
            'LineStyle','none','FaceColor',color.m3,'FaceAlpha',.1);
    view(3);
    xlabel('Time step size', 'fontweight', 'bold','Interpreter','Latex')
    ylabel('Samples', 'fontweight', 'bold','Interpreter','Latex')
    zlabel(sprintf('Comp Times for $X^{%d%d}$',row,col),...
        'fontweight', 'bold','Interpreter','Latex')
end
%% Output
    pic_type='eps';
    save_param='epsc';
    root=[pwd,'\','Pdf'];
    fileName=...
            sprintf('AB_const_T%1.3g_d%i_N%i_M%i',...
            param.T,param.d,param.N,param.M);
    tempPath=[root,'\','temp'];
    if exist(tempPath)==7
        [status, message, messageid] = rmdir(tempPath,'s');
    end
    mkdir(tempPath);
    inputPath=[root,'\','temp','\','input','.','tex'];
    delFile(inputPath);
    in_file=fopen(inputPath,'a');
    fprintf(in_file,...
        '\\section{Magnus expansion for $A$, $B$ constant and deterministic}\n');
    fprintf(in_file,...
        '\\subsection{Values of Moments}\n');
    for isample=1:1:length(sampleConv)
        for idelta=1:1:length(deltaConv)
            strMoment=sprintf('moment_Delta%g_M%d',...
                                      deltaConv(idelta),...
                                      sampleConv(isample));
            filePath=[root,'\','temp','\',strMoment,'.','tex'];
            fprintf(in_file,...
                        '\t\\input{%s}\n',changeSlash(filePath));
            file = fopen(filePath,'a');
            fprintf(file,...
                    '\\paragraph*{Configuration $\\Delta=%1.3g$, $M=%d$}\\hfill\\\\\n',...
                    deltaConv(idelta),...
                    sampleConv(isample));
            fprintf(file,...
                    '\\begin{tabular}{@{}*{6}{c}@{}}\n');
            fprintf(file,...
                    'Method & $E[(X^{11}_t)^k]$ & $E[(X^{12}_t)^k]$ & $E[(X^{21}_t)^k]$ & $E[(X^{22}_t)^k]$ & Total time\\\\\n');
            for mm=1:1:length(moments)
                morder=moments(mm);
                fprintf(file,'\\hline\n');
                fprintf(file,'\\multicolumn{6}{c}{$k=%d$}\\\\\n',morder);
                fprintf(file,'\\verb+euler+ & ');
                for row=1:1:param.d
                    for col=1:1:param.d
                        fprintf(file,'%g & ',momentEulerConv(row,col,mm,isample,idelta));
                    end
                end
                fprintf(file,'%g \\\\\n ',ctimeEulerConvTotal(mm,isample,idelta));
                fprintf(file,'\\verb+m1+ & ');
                for row=1:1:param.d
                    for col=1:1:param.d
                        fprintf(file,'%g & ',momentM1Conv(row,col,mm,isample,idelta));
                    end
                end
                fprintf(file,'%g \\\\\n ',ctimeM1ConvTotal(mm,isample,idelta));
                fprintf(file,'\\verb+m2+ & ');
                for row=1:1:param.d
                    for col=1:1:param.d
                        fprintf(file,'%g & ',momentM2Conv(row,col,mm,isample,idelta));
                    end
                end
                fprintf(file,'%g \\\\\n ',ctimeM2ConvTotal(mm,isample,idelta));
                fprintf(file,'\\verb+m3+ & ');
                for row=1:1:param.d
                    for col=1:1:param.d
                        fprintf(file,'%g & ',momentM3Conv(row,col,mm,isample,idelta));
                    end
                end
                if mm < length(moments)
                    fprintf(file,'%g \\\\\n',ctimeM3ConvTotal(mm,isample,idelta));
                else
                    fprintf(file,'%g \n',ctimeM3ConvTotal(mm,isample,idelta));
                end
            end
            fprintf(file,'\\end{tabular}\\hfill\\\\\n');
            fclose(file);
        end
    end
    fclose(in_file);
    templatePath=[pwd, '\' ,'template', '.','tex'];
    outputPath=[pwd,'\','template','.','pdf'];
%     tempPath=[root,'\','temp'];
    copyPath=[root,'\',fileName];
    copyPathFile=[root,'\',fileName,'\',fileName,'.','pdf'];
    if exist(copyPath)==7
        rmdir(copyPath,'s');
    end
    mkdir(copyPath);
%     delFile(copyPath);
    str1=sprintf('pdflatex %s',templatePath);
    system(str1);
    copyfile(tempPath,copyPath);
    copyfile(outputPath,copyPathFile);
%% Output
disp('done');
%% Aux functions
function str=changeSlash(str)
    for i=1:1:length(str)
        if strcmp(str(i),'\')
            str(i)='/';
        end
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function matrix2Latex(file,M)
    [n,m]=size(M);
    fprintf(file,...
        '\\left[\\begin{array}[c]{*{%i}{c}}\n',n);
    for i=1:1:n
        for j=1:1:m
            fprintf(file,...
                '%g ',M(i,j));
            if(j==m)
                fprintf(file,...
                    '\\\\\n');
            else
                fprintf(file,...
                    '&');
            end
        end
    end
    fprintf(file,...
        '\\end{array}\\right]\n');
end