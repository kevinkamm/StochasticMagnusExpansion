function Surfaces=AB_const_surface(Result,time,surfaces,varargin)
%% Plots predefined figures
%% Input:
% * see test_SPDE and SPDE_run
%% Output:
% * (1 x m) Figure
%%
% internal paramters
    fontsize=22;
    linewidth=1.5;
    facecolor.exact='b';
    facecolor.euler='r';
    facecolor.m1='b';
    facecolor.m2='y';
    facecolor.m3='c';
    show_title=1;
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'show_title'
                show_title=varargin{i+1};
        end
    end
    for kk=1:1:length(surfaces)
        Surfaces(kk)=figure('units','normalized',...
                    'outerposition',[0 0 1 1]); hold on;
        fig=surfaces{kk};
        figure_properties(Surfaces(kk));
        ij=fig{1};
        i=ij(1);
        j=ij(2);
        ind_t=fig{2};
        ind_p=fig{3};
        methods={fig{4:end}};
        [t,x]=...
            meshgrid(...
                reshape(time.t(ind_t),[length(ind_t) 1]),...
                reshape(ind_p,[length(ind_p) 1])...
            );
        for k=1:1:length(methods)
            Z=reshape(...
                Result.(methods{k}).X(i,j,ind_t,ind_p),...
                [length(ind_t),length(ind_p)]);
            s=surf(...
                x,...
                t,...
                Z',...
                'FaceAlpha',.5,...
                'FaceColor','interp',...
                'LineStyle','none'...
                );
            view(3);
            xlabel('paths', 'fontweight', 'bold')
            ylabel('t_i', 'fontweight', 'bold')
            str=sprintf('(X_t)^{%i,%i}',i,j);
            zlabel(str, 'fontweight', 'bold')
            if length(methods)==1
                c=parula(256);
                colormap(c);
            else
                set(s,'FaceColor',facecolor.(methods{k}));
            end
        end
        str=sprintf('X(%i,%i,:,:)',i,j);
        set(Surfaces(kk),'Name',str);
        if show_title
            title(str);
        end
        legend(methods,'Location','southoutside','NumColumns',length(methods));;
    end
    function figure_properties(fig)
        set(gca,'FontSize',fontsize)
        set(fig,'defaultlinelinewidth',linewidth)
        set(fig,'defaultaxeslinewidth',linewidth)
        set(fig,'defaultpatchlinewidth',linewidth)
        set(fig,'defaultAxesFontSize',fontsize)
    end
end