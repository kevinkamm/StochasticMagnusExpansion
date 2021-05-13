clear all;gpuDevice(1);close all

%% Parameters
% Gloabal
param.t0=0;
param.T=1;
param.d=2;
[A,B]=AB_const_coeff_initialize(param,'Example','default');
moments=[1];
% Z=2.576;% 99% confidence interval
Z=1.96; % 95% confidence interval
% Reference method
% %653.599 seconds
% device={'cpu','cpu'};
% deltaEulerRef=10^(4);
% sampleEulerRef = 10^5;
%272.599 seconds
device={'cpu','cpu'};
deltaEulerRef=100000+1;
sampleEulerRef = 10^4;
% device={'gpu','gpu'};
% deltaEulerRef=10000+1;
% sampleEulerRef = 10^4;
% Convergence test
deltaStart = 10+1;
deltaEnd = deltaEulerRef;
deltaConv = findN(deltaStart,deltaEnd);
deltaSteps = length(deltaConv);
sampleStart = sampleEulerRef;
sampleEnd = sampleEulerRef;
sampleSteps = 1;
sampleConv = linspace(sampleStart,sampleEnd,sampleSteps);
% Result variables
ctimeEulerRef=0;
ctimeEulerRefTotal=zeros(length(moments),1);
momentEulerRef=zeros(param.d,param.d,length(moments));
ctimeEulerConv=zeros(sampleSteps,deltaSteps);
ctimeEulerConvTotal=zeros(length(moments),sampleSteps,deltaSteps);
momentEulerConv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
rmseEulerConv=zeros(param.d,param.d,length(moments),sampleSteps,deltaSteps);
% Save paths
root=[pwd,'\','Euler','\','Matlab'];

%% Computing reference method
disp('Computing reference method')
% initialize
Nref=deltaEulerRef;
param.N=deltaEulerRef;
param.M=sampleEulerRef;
[time,BM]=AB_const_param_initialize(param);

% Calculate euler
disp('Calculate euler')
tic;
XeulerRef=AB_const_euler(time.dt,...
                    BM.dW,...
                    A,...
                    B,...
                    'Comp Device',device{1}, 'disp', 0);
ctimeEulerRef=toc;               
fprintf('Elapsed time %g\n',ctimeEulerRef);
% Calculate moments
disp('Calculate moments')
for mm=1:1:length(moments)
    morder=moments(mm);
    tic;
    momentEulerRef(:,:,mm)=mean(XeulerRef(:,:,end,:).^morder,4);
    ctimeEulerRefTotal(mm)=ctimeEulerRef+toc;
end
clear XeulerRef

%% Computing convergence test
for isample=1:1:length(sampleConv)
    for idelta=1:1:length(deltaConv)
        gpuDevice(1);
        fprintf('Runs convergence: %d/%d samples, %d/%d deltas\n',...
            isample,length(sampleConv),idelta,length(deltaConv))
        % initialize
        param.N=deltaConv(idelta);
        ind_t=[1:1:param.N];
        ind_t(2:1:end)=ind_t(1:1:end-1).*floor((Nref-1)/(param.N-1))+1;

        tref=linspace(param.t0,param.T,Nref);
        dtref=diff(tref);
        t=linspace(param.t0,param.T,param.N);
        fprintf('diff of t %g\n',sum(abs(reshape(time.t(1,1,ind_t,1),[],1)-reshape(t,[],1)),'all'))
        dt=(param.T-param.t0)/(param.N-1);
        dW=diff(BM.W(:,:,ind_t,:),1,3);
        % Calculate euler
        disp('Calculate euler')
        tic;
        X=AB_const_euler(time.dt,...
                            dW,...
                            A,...
                            B,...
                            'Comp Device',device{2});
        ctimeEulerConv(isample,idelta)=toc;
        fprintf('Elapsed time %g\n',ctimeEulerConv(isample,idelta));
        % Calculate moments for euler
        disp('Calculate moments for euler')
        for mm=1:1:length(moments)
            morder=moments(mm);
            tic;
            momentEulerConv(:,:,mm,isample,idelta)=mean(X(:,:,end,:).^morder,4);
            ctimeEulerConvTotal(mm,isample,idelta)=ctimeEulerConv(isample,idelta)+toc;
            % Calculate rmse for euler
            rmseEulerConv(:,:,mm,isample,idelta)=...
                sqrt(mean((X(:,:,end,:).^morder-reshape(momentEulerRef(:,:,mm),param.d,param.d,1,1)).^2,4)./size(X,4));
        end
    end
end
clear X
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
                plot(deltaConv,squeeze(rmseEulerConv(row,col,mm,isample,:)),...
                    'LineStyle','-','Color',color.euler);
            end
            xlabel('$\Delta$', 'fontweight', 'bold','Interpreter','Latex')
            ylabel(sprintf('RMSE for $X^{%d%d}$',row,col),...
                'fontweight', 'bold','Interpreter','Latex')
        end
    end
end
% moment
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
                plot(deltaConv,squeeze(abs(momentEulerConv(row,col,mm,isample,:)-...
                    momentEulerRef(row,col,mm)))',...
                    'LineStyle','-','Color',color.euler);
            end
            xlabel('$\Delta$', 'fontweight', 'bold','Interpreter','Latex')
            ylabel(sprintf('Abs Err for $X^{%d%d}$',row,col),...
                'fontweight', 'bold','Interpreter','Latex')
        end
    end
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
%% Output
disp('done');
%% Aux functions
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function N=findN(N0,N1)
    k=1;
    N=-1;
    for i=N0:1:N1
        if mod((N1-1),(i-1))==0
            N(k)=i;
            k=k+1;
        end
    end
end