function Plots=AB_const_plot(Result,time,plots,varargin)
%% Plots predefined figures
%% Input:
% * see test_AB_const and AB_const_run
%% Output:
% * (1 x m) Figure
%%
% internal paramters
    fontsize=22;
    linewidth=1;
%     linspecs.exact='k-';
    linspecs.exact.linestyle={'-','-'};
    linspecs.exact.marker={'none','none'};
    linspecs.exact.color={'k','k'};
%     linspecs.euler='r--';
    linspecs.euler.linestyle={'-.','-'};
    linspecs.euler.marker={'none','none'};
    linspecs.euler.color={'r',rgb2gray([1 0 0])};
%     linspecs.m1='b.';
    linspecs.m1.linestyle={'none','none'};
    linspecs.m1.marker={'.','.'};
    linspecs.m1.color={'b',rgb2gray([1 1 0])};
%     linspecs.m2='yo';
    linspecs.m2.linestyle={'none','none'};
    linspecs.m2.marker={'o','o'};
    linspecs.m2.color={'y',rgb2gray([0 1 0])};
%     linspecs.m3='gx';
    linspecs.m3.linestyle={'none','none'};
    linspecs.m3.marker={'x','x'};
    linspecs.m3.color={'g',rgb2gray([0 0 1])};
    show_title=1;
    resolution=size(time.t,3);
    blackwhite='both'; %false and true are possible too
    for kk=1:2:length(varargin)
        switch varargin{kk}
            case 'show_title'
                show_title=varargin{kk+1};
            case 'resolution'
                resolution=varargin{kk+1};
            case 'blackwhite'
                blackwhite=varargin{kk+1};
        end
    end
    if resolution >= size(time.t,3)
        res_step=1;
    else
        res_step=floor(size(time.t,3)/resolution);
    end
    t_pos=[1:res_step:size(time.t,3)];
    switch blackwhite
        case 'true'
            bw_plot(2);
        case 'both'
            bw_plot(1);
            bw_plot(2);
        otherwise
            bw_plot(1);
    end
    function bw_plot(bw)
        if exist('Plots') && ~isempty(Plots)
            ii=length(Plots);
        else
            ii=0;
        end
        for i=1:1:length(plots)
            Plots(ii+i)=figure('units','normalized',...
                        'outerposition',[0 0 1 1]); hold on;
            figure_properties(Plots(i));
            fig=plots{i};
            pos=fig{1};
            methods={fig{2:end}};
            for k=1:1:length(methods)
                plot(...
                    reshape(time.t(t_pos),[size(time.t(t_pos),3),1]),...
                    reshape(...
                        Result.(methods{k}).X(pos(1),pos(2),t_pos,pos(3)),...
                        [size(time.t(t_pos),3) 1]),...
                    'LineStyle',...
                    linspecs.(methods{k}).linestyle{bw},...
                    'Marker',...
                    linspecs.(methods{k}).marker{bw},...
                    'Color',...
                    linspecs.(methods{k}).color{bw}...
                    );
            end
            str=sprintf('(X_t)_{%i%i}',pos(1),pos(2));
            set(Plots(i),'Name',str);
            xlabel('t', 'fontweight', 'bold')
            ylabel(str, 'fontweight', 'bold')
            if show_title
                title(str);
            end
            legend(methods,'Location','southoutside','NumColumns',length(methods));  
        end
    end
    function figure_properties(fig)
        set(gca,'FontSize',fontsize)
        set(fig,'defaultlinelinewidth',linewidth)
        set(fig,'defaultaxeslinewidth',linewidth)
        set(fig,'defaultpatchlinewidth',linewidth)
        set(fig,'defaultAxesFontSize',fontsize)
    end
end