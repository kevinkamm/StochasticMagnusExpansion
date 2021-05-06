function [Result,varargout]=B0_var_run(time,BM,ind_t,ind_p,A,methods,varargin)
%% Evaluates for declared methods and given parameters the following SDE
% 
% $$ dX_t = A X_t dW_t, X_0=I_d $$
% 
%% Input:
% * see test_SPDE
%% Output:
% * (struct) Result: 
%%
% # (d x d x N x N array) Result.(method).X containing the result
% # (double) Result.(method).ctime.[total,expm,log] containing the computational time
%%
%
% if ~isempty(varargin)
%     for k=1:1:length(varargin)
%         switch varargin{k}
%         end
%     end
% end
%
for k=1:1:length(methods)
    method=methods{k};
    fprintf('Calculating %s\n',method{1});
    switch method{1}
        case 'exact'
            tic;
            X=...
                B0_var_exact(...
                    A.f11,...
                    A.f12,...
                    A.f22,...
                    time.t_fine,...
                    BM.W_fine,...
                    method{2}{:});
            Result.exact.ctime.total=toc;
            fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            Result.exact.X=X(:,:,ind_t,ind_p); clear X;
        case 'euler'
            temp=zeros(2,2,size(time.t_fine,3));
            temp(1,1,:)=A.f11;
            temp(1,2,:)=A.f12;
            temp(2,2,:)=A.f22;
            tic;
            X=...
                B0_var_euler(...
                    BM.dW_fine,...
                    temp,...
                    method{2}{:});
            Result.euler.ctime.total=toc;
            fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            Result.euler.X=X(:,:,ind_t,ind_p); clear X temp;
        case {'m1','m2','m3'}
            m=method{1};
            fprintf('Calculating log\n');
            tic;
            O=...
                B0_var_magnus(...
                    time.t,...
                    BM.W,...
                    A.f11(ind_t),...
                    A.f12(ind_t),...
                    A.f22(ind_t),...
                    A.df11(ind_t),...
                    A.df12(ind_t),...
                    A.df22(ind_t),...
                    str2num(m(end)),...
                    method{2}{:});
            Result.(method{1}).ctime.log=toc;
            fprintf('Elapsed time %g\n',Result.(method{1}).ctime.log);
            clear m;
            fprintf('Calculating matrix exp\n');
            tic;
            Result.(method{1}).X=m_exp(...
                O,...
                'type','matlab',...
                'Comp Device','cpuReshape');
            Result.(method{1}).ctime.expm=toc;
            fprintf('Elapsed time %g\n',Result.(method{1}).ctime.expm);
            clear O;
            Result.(method{1}).ctime.total=...
                Result.(method{1}).ctime.expm+...
                Result.(method{1}).ctime.log;
        otherwise
            disp('Unknown method');
    end
end
end