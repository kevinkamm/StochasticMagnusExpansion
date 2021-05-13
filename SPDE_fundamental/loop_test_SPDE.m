clear;gpuDevice(1);
for n=1:1:1
clc;fprintf('Loop %i\n',n);
%% Parameters
%% 
% # time axis: 
%%
% * $[\verb+t0+,\verb+T+]\subset \bf{R}_{\geq 0}$
% * homogeneous time grid with refinement $\verb+N+\in \bf{N}$: 
%       $\verb+t0+=t_{1}\leq \dots \leq t_{N}=\verb+T+$
% time scale: fix N_1 for T=1 to 1000, 10000
%   then: N_t= ceil(1000* t),ceil(10000 * t)
Delta=1/10000;
param.t0=0;
param.T=.5;
param.N=ceil(Delta^(-1)*param.T);
%%
% # sample paths of Brownian motion $\verb+M+\in\bf{N}$
param.M=50;
%%
% # spatial domain:
%%
% * cut-off: $[\verb+xa+,\verb+xb+]\subset \bf{R}$
% * discretization of spatial domain with $d+2$ points: 
%           $xa=x_1\leq\dots\leq x_{d+2}=xb$
param.d=51;
param.xa=-2;
param.xb=2;
param.kappa=floor(param.d/2);
%% SPDE
%%
% 
% $$ du(t,x)=\frac{a}{2} (\partial_{xx} u)(t,x) dt + \sigma (\partial_{x} u)(t,x) dW_t $$
% with a> \sigma^2
% 
SPDE.sigma=.15;
SPDE.a=.2;
%% Initialize
[time,space,BM]=SPDE_param_initialize(param);
SPDE=SPDE_coeff_initialize(param,space,SPDE);
% significance level:
i01=floor((param.d-param.kappa)./2);
i02=i01+param.kappa;
%% Method Selection
% methods={method_1,...,method_n}
% one method consists of: method_k={name,{varargin}}
% available method names: exact, euler, m1, m2, m3
methods=...
    {...
        {'exact',{'Comp Device','cpu'}},...
        {'euler',{'Comp Device','cpu'}},...
        {'m1',{'Comp Device','cpu'}},...
        {'m2',{'Comp Device','cpu'}},...
        {'m3',{'Comp Device','cpu'}}...
    };
%% Plot Selection
% two-dimensional plots:
% solution path plots:
% plots={plot_1,...,plot_n}
% one plot consists of: plot_k={[i j m],method_1,...,method_j}, where
%   m is the path of the Brownian motion
%   (i,j) corresponds to an indices in the solution matrix X_t(i,j)
%   e.g. i=j are the diagonal entries
% plots=...
%     {...
%         {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m1','m2','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m2','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m1','m3'}...
%     };
% error plots (ecdf):
% error_plots={err_plot_1,...,err_plot_n}
% one plot consists of: err_plot_k={sp_shape,method_r,{method_1,...,method_j}}, where
%   sp_shape is the subplot shape, e.g. sp_shape=[2 3] for 2 rows and 
%       3 columns of subplots
%   method_r is the reference method for the errors, e.g. 'euler' 
% error_plots=...
%     {...
%         {[2 2],'exact',{'euler'}},...
%         {[2 2],'exact',{'m1'}},...
%         {[2 2],'exact',{'m2'}},...
%         {[2 2],'exact',{'m3'}}...
%     };
% three-dimensional plots:
% surfaces={surface_1,...,surface_n}
% one surface consists of: surface_k={[i,j],[s],[m],method_1,...,method_z}
%   (i,j) corresponds to the indizes in the solution matrix X_t(i,j)
%   m corresponds to the indizes in the paths
%   s corresponds to the indizes in the time grid
% surfaces=...
%     {...
%         {[floor(param.d/2) floor(param.d/2)],[1:1:param.N],[1:1:param.M],'euler'},...
%         {[floor(param.d/2) floor(param.d/2)],[1:1:param.N],[1:1:param.M],'m3'},...
%         {[floor(param.d/2) floor(param.d/2)],[1:1:param.N],[1:1:param.M],'euler','m3'},...
%         {[floor(param.d/2) floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'euler'},...
%         {[floor(param.d/2) floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'m3'},...
%         {[floor(param.d/2) floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'euler','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)],[1:1:param.N],[1:1:param.M],'euler'},...
%         {[floor(param.d/2)+1 floor(param.d/2)],[1:1:param.N],[1:1:param.M],'m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)],[1:1:param.N],[1:1:param.M],'euler','m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'euler'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'m3'},...
%         {[floor(param.d/2)+1 floor(param.d/2)+1],[1:1:param.N],[1:1:param.M],'euler','m3'}...
%     };
%% Error Selection
% path-vise errors
% errors={errror_1,...,error_n}
% one error consists of: error_k={i,j,method_r,{method_1,...,method_u}}, where
%   (i,j) corresponds to an indices in the solution matrix X_t(i,j)
%   and method_r is the reference method
% errors=...
%     {...
%         {floor(param.d/2),floor(param.d/2),'exact',{'euler','m1','m2','m3'}}...
%     };
% total errors
% err_total={err_total_1,...,err_total_2}
% one error consists of: err_total_k={method_r,{method_1,...,method_u}}
errors_total=...
    {...
        {'exact',{'euler','m1','m2','m3'}}...
    };
%% Output Selection
% output={output_1,...,output_n}
% one output consists of output_k={name,varargin}
% available output: pdf
output=...
    {...
        {'Matlab',n}...
    };
%% Calculate methods
Result=SPDE_run(param,time,space,BM,SPDE,methods,'memory',0);
%% Calculate errors
if exist('errors')==1 && ~isempty(errors)
    Errors=SPDE_errors(Result,time,errors,'p',1,'T',param.T);
end
if exist('errors_total')==1 && ~isempty(errors_total)
    Errors.total=SPDE_errors_total(Result,param,time,[i01:1:i02],errors_total,'percentile',.99);
end
%% Plot 2D Figures
if exist('plots')==1 && ~isempty(plots)
    Plots=SPDE_plot(Result,time,plots,'show_title',0);
end
if exist('error_plots')==1 && ~isempty(error_plots)
    Error_Plots=SPDE_plot_err(Errors,error_plots,'show_title',0,'percentile',.99);
end
%% Plot 3D Figures
if exist('surfaces')==1 && ~isempty(surfaces)
    Surfaces=SPDE_surface(Result,time,space,surfaces,'show_title',0);
end
%% Create output
if exist('output')==1 && ~isempty(output)
    vararg={};
    if exist('Errors')==1 && ~isempty(Errors)
        vararg{end+1}='Errors';
        vararg{end+1}=Errors;
    end
    if exist('Plots')==1 && ~isempty(Plots)
        vararg{end+1}='Plots';
        vararg{end+1}=Plots;
    end
    if exist('Error_Plots')==1 && ~isempty(Error_Plots)
        vararg{end+1}='Error_Plots';
        vararg{end+1}=Error_Plots;
    end
    if exist('Surfaces')==1 && ~isempty(Surfaces)
        vararg{end+1}='Surfaces';
        vararg{end+1}=Surfaces;
    end
    SPDE_output(param,SPDE,Result,output,vararg{:});
end
end
%% Plots
clc;
samples=n*param.M;
Err.euler.E=zeros(1,samples);
Err.m1.E=zeros(1,samples);
Err.m2.E=zeros(1,samples);
Err.m3.E=zeros(1,samples);
root=[pwd,'\','Matlab'];
fileName_old=...
        sprintf('SPDE_T%1.3g_d%i_N%i_M%i_xa%i_xb_%i',...
        param.T,param.d,param.N,param.M,param.xa,param.xb);
fileName=...
        sprintf('SPDE_T%1.3g_d%i_N%i_M%i_xa%i_xb_%i',...
        param.T,param.d,param.N,samples,param.xa,param.xb);
if exist([root,'\',fileName_old])==7 && ~strcmp(fileName_old,fileName)
    mkdir([root,'\',fileName]);
    movefile([root,'\',fileName_old,'\','*'],[root,'\',fileName]);
    rmdir([root,'\',fileName_old]);
end
for m=1:1:n
    file=sprintf('Errors_%i.mat',m);
    S=load([root,'\',fileName,'\',file]);
    references=fieldnames(S);
    for ii=1:1:length(references)
        methods=fieldnames(S.(references{ii}));
        for jj=1:1:length(methods)
            Err_t=S.(references{ii}).(methods{jj}).Err;
            for kk=1:1:size(Err_t,2)
                E=Err.(methods{jj}).E;
                E(1,(m-1)*param.M+kk)=Err_t(kk);
                Err.(methods{jj}).E=E;
            end
        end
    end
end
fontsize=22;
linewidth=1.5;
%     linspecs.exact='k-';
linspecs.exact.linestyle={'-','-'};
linspecs.exact.marker={'none','none'};
linspecs.exact.color={'k','k'};
%     linspecs.euler='r--';
linspecs.euler.linestyle={'-.','none'};
linspecs.euler.marker={'none','+'};
linspecs.euler.color={'r',rgb2gray([1 0 0])};
%     linspecs.m1='b.';
linspecs.m1.linestyle={'none','none'};
linspecs.m1.marker={'.','.'};
linspecs.m1.color={'b',rgb2gray([0 0 1])};
%     linspecs.m2='yo';
linspecs.m2.linestyle={'none','none'};
linspecs.m2.marker={'o','o'};
linspecs.m2.color={'y',rgb2gray([1 1 0])};
%     linspecs.m3='gx';
linspecs.m3.linestyle={'none','none'};
linspecs.m3.marker={'x','x'};
linspecs.m3.color={'g',rgb2gray([0 1 0])};
percentile=.99;
bw=1;
pic_type='eps';
save_param='epsc';
pplots=fieldnames(Err);
for m=1:1:length(pplots)
    fig(m)=figure('units','normalized',...
                        'outerposition',[0 0 1 1]);
    set(gca,'FontSize',fontsize)
    set(fig(m),'defaultlinelinewidth',linewidth)
    set(fig(m),'defaultaxeslinewidth',linewidth)
    set(fig(m),'defaultpatchlinewidth',linewidth)
    set(fig(m),'defaultAxesFontSize',fontsize)
    subplot(2,2,1)
    Err.(pplots{m}).Exp=mean(Err.(pplots{m}).E);
    fprintf('%s:\n\tErr_t: %g\n',pplots{m},mean(Err.(pplots{m}).E));
    [f,x]=ecdf(Err.(pplots{m}).E);
    c=find(f>percentile,1,'first');
    if ~isempty(c)
        f=f(1:1:c);
        x=x(1:1:c);
    end
    plot(x,f,...
        'LineStyle',...
        linspecs.(pplots{m}).linestyle{bw},...
        'Marker',...
        linspecs.(pplots{m}).marker{bw},...
        'Color',...
        linspecs.(pplots{m}).color{bw});
    legend(pplots{m},'Location','southeast');
    str=sprintf('CDF_{Err_{%2.3g}}(x)',...
        param.T);
%     ytickformat('%1.2f')
    a=[cellstr(num2str(get(gca,'xtick')'*100))]; 
    pct = char(ones(size(a,1),1)*'%');
    new_xticks = [char(a),pct];
    set(gca,'xticklabel',new_xticks)
    a=[cellstr(num2str(get(gca,'ytick')'))]; 
    a{end}=char(num2str(percentile));
    new_yticks = [char(a)];
    set(gca,'yticklabel',new_yticks)
    xlabel('x', 'fontweight', 'bold')
    ylabel(str, 'fontweight', 'bold')
    fig_name=...
        sprintf('error_plot_%i',m);
    pngPath = [root,'\',fileName,'\',fig_name,'.',pic_type];
    if exist(pngPath)
        delete(pngPath);
    end
    saveas(fig(m),pngPath,save_param);
%     pngPath = [root,'\',fileName,'\',fig_name,'.','pdf'];
%     if exist(pngPath)
%         delete(pngPath);
%     end
%     saveas(fig(m),pngPath);
end
mean_name=...
        sprintf('Errors');
meanPath = [root,'\',fileName,'\',mean_name,'.','mat'];
if exist(meanPath)
    delete(meanPath);
end
save(meanPath,'-struct','Err');