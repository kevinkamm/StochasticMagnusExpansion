function [Result,varargout]=AB_const_run(time,BM,ind_t,ind_p,A,B,methods,varargin)
%% Evaluates for declared methods and given parameters the following SDE
% 
% $$ dX_t = B X_t dt + A X_t dW_t, X_0=I_d $$
% 
%% Input:
% * see test_SPDE
%% Output:
% * (struct) Result: 
%%
% # ($d\times N \times M$ array) Result.(method).u containing the result
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
%%
% # memory profiling
% memProfiling = 1;
for k=1:1:length(methods)
    method=methods{k};
    fprintf('Calculating %s\n',method{1});
%     if memProfiling
%         profile -memory on
%         profile clear
%     end
    switch method{1}
%         case 'exact'
%             tic;
%             Result.exact.u=...
%                 SPDE_exact_u0(...
%                     SPDE.a,...
%                     SPDE.sigma,...
%                     time.t,...
%                     space.x,...
%                     BM.W,...
%                     u0.u,...
%                     param.K,...
%                     method{2}{:});
%             Result.exact.ctime.total=toc;
%             fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
        case 'euler'
            tic;
            X=...
                AB_const_euler(...
                    time.dt_fine,...
                    BM.dW_fine,...
                    A,...
                    B,...
                    method{2}{:});
            Result.euler.ctime.total=toc;
            fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            Result.euler.X=X(:,:,ind_t,ind_p); clear X;
        case {'m1','m2','m3'}
            m=method{1};
            fprintf('Calculating log\n');
            tic;
            O=...
                AB_const_magnus(...
                    time.t,...
                    BM.W,...
                    A,...
                    B,...
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
%     if memProfiling
%         profile off
%         profsave(profile('info'),sprintf('memProf_%s',method{1}))
%     end
end
end