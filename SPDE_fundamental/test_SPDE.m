clear;gpuDevice(1);
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
param.T=.1;
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
param.d=101;
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
plots=...
    {...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m1','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'exact','euler','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m1','m3'}...
    };
% error plots (ecdf):
% error_plots={err_plot_1,...,err_plot_n}
% one plot consists of: err_plot_k={sp_shape,method_r,{method_1,...,method_j}}, where
%   sp_shape is the subplot shape, e.g. sp_shape=[2 3] for 2 rows and 
%       3 columns of subplots
%   method_r is the reference method for the errors, e.g. 'euler' 
error_plots=...
    {...
        {[2 2],'exact',{'euler'}},...
        {[2 2],'exact',{'m1'}},...
        {[2 2],'exact',{'m2'}},...
        {[2 2],'exact',{'m3'}}...
    };
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
        {'pdf'}...
    };
%% Calculate methods
Result=SPDE_run(param,time,space,BM,SPDE,methods,'memory',20);
% Result=SPDE_run(param,time,space,BM,SPDE,methods);
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
disp('done');