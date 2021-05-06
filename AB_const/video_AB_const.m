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
Delta=1/1000;
Delta_fine=1/10000;
param.t0=0;
param.T=.75;
param.N_fine=ceil(Delta_fine^(-1)*param.T);
param.N=ceil(Delta^(-1)*param.T);
%%
% # sample paths of Brownian motion $\verb+M+\in\bf{N}$
param.M_fine=10000;
param.M=10000;
%%
% # spatial domain:
% %
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
[time,BM,ind_t,ind_p]=AB_const_param_initialize(param);
[A,B]=AB_const_coeff_initialize(param,'Example','default');
%% Calculate methods
Result=AB_const_run(time,BM,ind_t,ind_p,A,B,methods);
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
    Plots=AB_const_plot(Result,time,plots,'show_title',0);
end
if exist('error_plots')==1 && ~isempty(error_plots)
    Error_Plots=AB_const_plot_err(Errors,error_plots,'show_title',0,'percentile',.99);
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
%% Videos
EulerColor={{'Color','r','LineStyle','none','Marker','*'},...
    {'Color','r','LineStyle','-','Marker','none'}};
M1Color={{'Color','y','LineStyle','none','Marker','*'},...
    {'Color','y','LineStyle','--','Marker','none'}};
M2Color={{'Color','g','LineStyle','none','Marker','*'},...
    {'Color','g','LineStyle','--','Marker','none'}};
M3Color={{'Color','b','LineStyle','none','Marker','*'},...
    {'Color','b','LineStyle','--','Marker','none'}};
for ii=1:1:param.d
    for jj=1:1:param.d
        close all;
% Video
fontsize=22;
linewidth=1.5;
t=reshape(time.t,[param.N 1]);
Euler=reshape(Result.euler.X(ii,jj,:,1),[param.N 1]);
M1=reshape(Result.m1.X(ii,jj,:,1),[param.N 1]);
M2=reshape(Result.m2.X(ii,jj,:,1),[param.N 1]);
M3=reshape(Result.m3.X(ii,jj,:,1),[param.N 1]);

%%first half
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)

t1=find(t>0.3,1,'first');ind1=1:1:t1-1;
minEuler=min(Euler(ind1));maxEuler=max(Euler(ind1));
minM1=min(M1(ind1));maxM1=max(M1(ind1));
minM2=min(M2(ind1));maxM2=max(M2(ind1));
minM3=min(M3(ind1));maxM3=max(M3(ind1));
minTotal=min([minEuler,minM1,minM2,minM3]);
maxTotal=max([maxEuler,maxM1,maxM2,maxM3]);
axis([t(1),t(t1-1),minTotal,maxTotal]);
subfolder=sprintf('X%i%i',ii,jj);
videopath=['Video','/',subfolder];
% Video 1
v=VideoWriter([videopath,'/','plot1'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),Euler(i),EulerColor{1}{:});hold on;
    legend('euler','Location','southoutside','NumColumns',1)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(t1-1),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
legend('euler','Location','southoutside','NumColumns',1)
saveas(fig,[videopath,'/','plot1.jpg']);
% Video 2
v=VideoWriter([videopath,'/','plot2'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M1(i),M1Color{1}{:});hold on;
    legend('euler','m1','Location','southoutside','NumColumns',2)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(t1-1),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
legend('euler','m1','Location','southoutside','NumColumns',2)
saveas(fig,[videopath,'/','plot2.jpg']);
% Video 3
v=VideoWriter([videopath,'/','plot3'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M2(i),M2Color{1}{:});hold on;
    legend('euler','m1','m2','Location','southoutside','NumColumns',3)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(t1-1),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
plot(t(ind1),M2(ind1),M2Color{2}{:});hold on;
legend('euler','m1','m2','Location','southoutside','NumColumns',3)
saveas(fig,[videopath,'/','plot3.jpg']);
% Video 4
v=VideoWriter([videopath,'/','plot4'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M3(i),M3Color{1}{:});hold on;
    legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(t1-1),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
plot(t(ind1),M2(ind1),M2Color{2}{:});hold on;
plot(t(ind1),M3(ind1),M3Color{2}{:});hold on;
legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
saveas(fig,[videopath,'/','plot4.jpg']);
close all;
% second half
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)

t1=find(t>0.3,1,'first');ind1=t1:1:param.N;
minEuler=min(Euler(ind1));maxEuler=max(Euler(ind1));
minM1=min(M1(ind1));maxM1=max(M1(ind1));
minM2=min(M2(ind1));maxM2=max(M2(ind1));
minM3=min(M3(ind1));maxM3=max(M3(ind1));
minTotal=min([minEuler,minM1,minM2,minM3]);
maxTotal=max([maxEuler,maxM1,maxM2,maxM3]);
axis([t(t1),t(end),minTotal,maxTotal]);

% Video 1
v=VideoWriter([videopath,'/','plot5'],'MPEG-4');
v.FrameRate=40;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),Euler(i),EulerColor{1}{:});hold on;
    legend('euler','Location','southoutside','NumColumns',1)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(t1),t(end),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
legend('euler','Location','southoutside','NumColumns',1)
saveas(fig,[videopath,'/','plot5.jpg']);
% Video 2
v=VideoWriter([videopath,'/','plot6'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M1(i),M1Color{1}{:});hold on;
    legend('euler','m1','Location','southoutside','NumColumns',2)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(t1),t(end),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
legend('euler','m1','Location','southoutside','NumColumns',2)
saveas(fig,[videopath,'/','plot6.jpg']);
% Video 3
v=VideoWriter([videopath,'/','plot7'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M2(i),M2Color{1}{:});hold on;
    legend('euler','m1','m2','Location','southoutside','NumColumns',3)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(t1),t(end),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
plot(t(ind1),M2(ind1),M2Color{2}{:});hold on;
legend('euler','m1','m2','Location','southoutside','NumColumns',3)
saveas(fig,[videopath,'/','plot7.jpg']);
% Video 4
v=VideoWriter([videopath,'/','plot8'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M3(i),M3Color{1}{:});hold on;
    legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(t1),t(end),minTotal,maxTotal]);
plot(t(ind1),Euler(ind1),EulerColor{2}{:});hold on;
plot(t(ind1),M1(ind1),M1Color{2}{:});hold on;
plot(t(ind1),M2(ind1),M2Color{2}{:});hold on;
plot(t(ind1),M3(ind1),M3Color{2}{:});hold on;
legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
saveas(fig,[videopath,'/','plot8.jpg']);
close all;
% Total 
ind1=1:1:param.N;
minEuler=min(Euler(ind1));maxEuler=max(Euler(ind1));
minM1=min(M1(ind1));maxM1=max(M1(ind1));
minM2=min(M2(ind1));maxM2=max(M2(ind1));
minM3=min(M3(ind1));maxM3=max(M3(ind1));
minTotal=min([minEuler,minM1,minM2,minM3]);
maxTotal=max([maxEuler,maxM1,maxM2,maxM3]);
%Video total 1
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(end),minTotal,maxTotal]);
v=VideoWriter([videopath,'/','total1'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),Euler(i),EulerColor{1}{:});hold on;
    legend('euler','Location','southoutside','NumColumns',1)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(end),minTotal,maxTotal]);
plot(t,Euler,EulerColor{2}{:});hold on;
legend('euler','Location','southoutside','NumColumns',1)
saveas(fig,[videopath,'/','total1.jpg']);
%Video total 2
v=VideoWriter([videopath,'/','total2'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M1(i),M1Color{1}{:});hold on;
    legend('euler','m1','Location','southoutside','NumColumns',2)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(end),minTotal,maxTotal]);
plot(t,Euler,EulerColor{2}{:});hold on;
plot(t,M1,M1Color{2}{:});hold on;
legend('euler','m1','Location','southoutside','NumColumns',2)
saveas(fig,[videopath,'/','total2.jpg']);
%Video total 3
v=VideoWriter([videopath,'/','total3'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M2(i),M2Color{1}{:});hold on;
    legend('euler','m1','m2','Location','southoutside','NumColumns',3)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(end),minTotal,maxTotal]);
plot(t,Euler,EulerColor{2}{:});hold on;
plot(t,M1,M1Color{2}{:});hold on;
plot(t,M2,M2Color{2}{:});hold on;
legend('euler','m1','m2','Location','southoutside','NumColumns',3)
saveas(fig,[videopath,'/','total3.jpg']);
%Video total 4
v=VideoWriter([videopath,'/','total4'],'MPEG-4');
v.FrameRate=20;
v.Quality=100;
open(v)
for i=ind1
    plot(t(i),M3(i),M3Color{1}{:});hold on;
    legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
    drawnow;
    F=getframe(gcf);
    writeVideo(v,F);
end
close(v);close all;
fig=figure('units','normalized',...
           'outerposition',[0 0 1 1]); hold on;
set(gca,'FontSize',fontsize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
axis([t(1),t(end),minTotal,maxTotal]);
plot(t,Euler,EulerColor{2}{:});hold on;
plot(t,M1,M1Color{2}{:});hold on;
plot(t,M2,M2Color{2}{:});hold on;
plot(t,M3,M3Color{2}{:});hold on;
legend('euler','m1','m2','m3','Location','southoutside','NumColumns',4)
saveas(fig,[videopath,'/','total4.jpg']);
    end
end