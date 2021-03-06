clear;close all 
gpuDevice(1);
% rng('default');
% profile -memory on;
%% Parameters
%%
% loops for computational time
loop=1;
%% 
% # time axis: 
%%
% * $[\verb+t0+,\verb+T+]\subset \bf{R}_{\geq 0}$
% * homogeneous time grid with refinement $\verb+N+\in \bf{N}$: 
%       $\verb+t0+=t_{1}\leq \dots \leq t_{N}=\verb+T+$
% time scale: fix N_1 for T=1 to 1000, 10000
%   then: N_t= ceil(1000* t),ceil(10000 * t)
% Delta_fine=10^(-3);
% Delta=sqrt(Delta_fine);
Delta_fine=10^(-4);
Delta=10^(-2);
param.t0=0;
param.T=1;
%(param.N_fine-1)/(param.N-1)\in \mathbb{N}! 
param.N_fine=ceil(Delta_fine^(-1)*param.T)+1;
param.N=ceil(Delta^(-1)*param.T)+1;
if mod((param.N_fine-1),(param.N-1))~=0
    fprintf('Matching grids are given if (param.N_fine-1)/(param.N-1)\\in \\mathbb{N}\n')
    fprintf('Current mod((param.N_fine-1),(param.N-1))=%g\n',...
        mod((param.N_fine-1),(param.N-1)))
end
%%
% # sample paths of Brownian motion $\verb+M+\in\bf{N}$
param.M_fine=10^3;
param.M=param.M_fine; % Don't change
%%
% # spatial domain:
%%
% * dimension of matrices
param.d=2;
%% Method Selection
% methods={method_1,...,method_n}
% one method consists of: method_k={name,{varargin}}
% available method names: euler, m1, m2, m3
methods=...
    {...
        {'euler',{'Comp Device','gpu'}},...
        {'m1',{'Comp Device','gpu'}},...
        {'m2',{'Comp Device','gpu'}},...
        {'m3',{'Comp Device','gpu'}}...
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
        {[floor(param.d/2) floor(param.d/2) 1],'euler','m1','m2','m3'},...
        {[floor(param.d/2) floor(param.d/2) 1],'euler','m2','m3'},...
        {[floor(param.d/2) floor(param.d/2) 1],'euler','m3'},...
        {[floor(param.d/2) floor(param.d/2)+1 1],'euler','m1','m2','m3'},...
        {[floor(param.d/2) floor(param.d/2)+1 1],'euler','m2','m3'},...
        {[floor(param.d/2) floor(param.d/2)+1 1],'euler','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2) 1],'euler','m1','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2) 1],'euler','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2) 1],'euler','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m1','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m2','m3'},...
        {[floor(param.d/2)+1 floor(param.d/2)+1 1],'euler','m3'}...
    };
% error plots (ecdf):
% error_plots={err_plot_1,...,err_plot_n}
% one plot consists of: err_plot_k={sp_shape,method_r,{method_1,...,method_j}}, where
%   sp_shape is the subplot shape, e.g. sp_shape=[2 3] for 2 rows and 
%       3 columns of subplots
%   method_r is the reference method for the errors, e.g. 'euler' 
error_plots=...
    {...
        {[2 2],'euler',{'m1'}},...
        {[2 2],'euler',{'m2'}},...
        {[2 2],'euler',{'m3'}},...
        {[1 3],'euler',{'m1','m2','m3'}},...
        {[3 1],'euler',{'m1','m2','m3'}}...
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
errors=...
    {...
        {floor(param.d/2),floor(param.d/2),'euler',{'m1','m2','m3'}},...
        {floor(param.d/2),floor(param.d/2)+1,'euler',{'m1','m2','m3'}},...
        {floor(param.d/2)+1,floor(param.d/2),'euler',{'m1','m2','m3'}},...
        {floor(param.d/2)+1,floor(param.d/2)+1,'euler',{'m1','m2','m3'}}...
    };
% total errors
% err_total={err_total_1,...,err_total_2}
% one error consists of: err_total_k={method_r,{method_1,...,method_u}}
errors_total=...
    {...
        {'euler',{'m1','m2','m3'}}...
    };
%% Output Selection
% output={output_1,...,output_n}
% one output consists of output_k=name
% available output: pdf
output=...
    {...
        'pdf'...
    };
%% Initialize
tic;
[time,BM,ind_t,ind_p]=AB_const_param_initialize(param);
[A,B]=AB_const_coeff_initialize(param,'Example','default');
toc;
%% Calculate methods
Result=AB_const_run(time,BM,ind_t,ind_p,A,B,methods,'loop',loop);
%% Calculate errors
if exist('errors')==1 && ~isempty(errors)
%     Errors=AB_const_errors(Result,time,errors,'p',1,'T',param.T);
    Errors=AB_const_errors(Result,time,errors);
end
if exist('errors_total')==1 && ~isempty(errors_total)
    Errors.total=AB_const_errors_total(Result,param,time,errors_total,'percentile',.99);
end
%% Plot 2D Figures
if exist('plots')==1 && ~isempty(plots)
    Plots=AB_const_plot(Result,time,plots,'show_title',0,'blackwhite','false');
end
if exist('error_plots')==1 && ~isempty(error_plots)
    Error_Plots=AB_const_plot_err(Errors,error_plots,'show_title',0,'percentile',.99,'blackwhite','false');
end
%% Plot 3D Figures
if exist('surfaces')==1 && ~isempty(surfaces)
    Surfaces=AB_const_surface(Result,time,surfaces,'show_title',0);
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
    AB_const_output(param,A,B,Result,output,vararg{:});
end
disp('done');