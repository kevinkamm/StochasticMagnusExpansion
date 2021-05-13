function Error_Plots=AB_const_plot_err(Errors,error_plots,varargin)
    show_title=1;
    percentile=1;
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
    linspecs.m1.color={'b',rgb2gray([1 1 0])};
%     linspecs.m2='yo';
    linspecs.m2.linestyle={'none','none'};
    linspecs.m2.marker={'o','o'};
    linspecs.m2.color={'y',rgb2gray([0 1 0])};
%     linspecs.m3='gx';
    linspecs.m3.linestyle={'none','none'};
    linspecs.m3.marker={'x','x'};
    linspecs.m3.color={'g',rgb2gray([0 0 1])};
    blackwhite='both'; %false and true are possible too
    for kk=1:2:length(varargin)
        switch varargin{kk}
            case 'show_title'
                show_title=varargin{kk+1};
            case 'percentile'
                percentile=varargin{kk+1};
            case 'blackwhite'
                blackwhite=varargin{kk+1};
        end
    end
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
        if exist('Error_Plots') && ~isempty(Error_Plots)
            ii=length(Error_Plots);
        else
            ii=0;
        end
        for i=1:1:length(error_plots)
            Error_Plots(ii+i)=figure('units','normalized',...
                        'outerposition',[0 0 1 1]); hold on;
            figure_properties(Error_Plots(i));
            fig=error_plots{i};
            sp_shape=fig{1};
            ref=fig{2};
            methods=fig{3};
            for j=1:1:length(methods)
                subplot(sp_shape(1),sp_shape(2),j)
                f=Errors.total.(ref).(methods{j}).ecdf.f;
                x=Errors.total.(ref).(methods{j}).ecdf.x;
                plot(x,f,...
                    'LineStyle',...
                    linspecs.(methods{j}).linestyle{bw},...
                    'Marker',...
                    linspecs.(methods{j}).marker{bw},...
                    'Color',...
                    linspecs.(methods{j}).color{bw});
                legend(methods{j},'Location','southeast');
                str=sprintf('CDF_{Err_{%2.3g}}(x)',...
                    Errors.total.(ref).(methods{j}).ecdf.T);
                ytickformat('%1.2f')
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
            end
            str=sprintf('ecdf');
            set(Error_Plots(i),'Name',str);
            if show_title
                title(str);
            end
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