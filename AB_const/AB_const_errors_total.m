function [Total,varargout]=AB_const_errors_total(Result,param,time,errors_total,varargin)
%% Evaluates Errors of types: (absolute) L^p and (relative) L^1 for given Results
%% Input:
% * see test_AB_const and AB_const_run
%% Output:
% * (struct) Errors: Errors.(ip).(reference).(method).*
%%
% # ip: ('i%d_%d') where %d corresponds e.g. to the path i=151, j=51
% # reference: ('name') of reference method for error analysis
% # method: ('name') of method compared to reference for error analysis
% # * : either ('l%d_error') with default %d=2 or ('l1_error_rel')
% # (double) ('l%d_error') containing the L^p([0,T]\times \Omega,\bf(R)) error
% # (struct) ('l1_error_rel') containing fields: min, max, mean
%%
% internal paramters
% fontsize=22;
% linewidth=1.5;
norm_type='fro'; 
cutlvl=0; % cut off between zero and cutlvl for relative errors
percentile=100000;
%%
%
if ~isempty(varargin)
    for k=1:1:length(varargin)
        switch varargin{k}
            case 'norm_type'
                norm_type=varargin{k+1};
            case 'cutlvl'
                cutlvl=varargin{k+1};
            case 'percentile'
                percentile=varargin{k+1};
        end
    end
end
varargout{1}=norm_type;
varargout{2}=cutlvl;
varargout{3}=percentile;
for k=1:1:length(errors_total)
    error=errors_total{k};
    fprintf('Displaying Total Errors:\n');
    reference=error{1};
    methods=error{2};
    X_ref=Result.(reference).X;
    fprintf('\t Reference method: %s\n',reference);
    for i=1:1:length(methods)
        fprintf('\t\t %s:\n',methods{i});
        X=Result.(methods{i}).X;
        %Error time and sample dependent
        Err_t=err_t(X_ref,X);
        %Expectation at terminal time
        Total.(reference).(methods{i}).expectation=...
            mean(Err_t(end,:));
        fprintf('\t\t\t E(err_T):\n\t\t\t\t %g\n',...
            Total.(reference).(methods{i}).expectation);
        %Empirical CDF at terminal time
        [f,x]=ecdf(Err_t(end,:));
        c=find(f>percentile,1,'first');
        if ~isempty(c)
            f=f(1:1:c);
            x=x(1:1:c);
        end
        Total.(reference).(methods{i}).ecdf.f=f;
        Total.(reference).(methods{i}).ecdf.x=x;
        Total.(reference).(methods{i}).ecdf.T=param.T;
%             %plotting cdf
%             CDFPLOT=...
%                 figure('units','normalized',...
%                         'outerposition',[0 0 1 1]);
%             figure_properties(CDFPLOT);
%             plot(x,f);
%             Total.(reference).(methods{i}).CDFPLOT=CDFPLOT;
%             %histogram
%             ECDFhist=...
%                 figure('units','normalized',...
%                         'outerposition',[0 0 1 1]);
%             figure_properties(ECDFhist);
%             ecdfhist(f,x);
%             Total.(reference).(methods{i}).ECDFhist=ECDFhist;
    end
end
    function Err_t=err_t(X_ref,X)
        Err_t=zeros(size(X_ref,3),size(X_ref,4));
        for ii=1:1:size(X_ref,3)
            for jj=1:1:size(X_ref,4)
                Err_t(ii,jj)=...
                    (norm(X_ref(:,:,ii,jj)-X(:,:,ii,jj),norm_type)./...
                    norm(X_ref(:,:,ii,jj),norm_type));
            end
        end
        Err_t=cumsum(Err_t.*time.dt,1)./reshape(time.t,[size(X_ref,3) 1]);
    end
%     function figure_properties(fig)
%         set(gca,'FontSize',fontsize)
%         set(fig,'defaultlinelinewidth',linewidth)
%         set(fig,'defaultaxeslinewidth',linewidth)
%         set(fig,'defaultpatchlinewidth',linewidth)
%         set(fig,'defaultAxesFontSize',fontsize)
%     end
end