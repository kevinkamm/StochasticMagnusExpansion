function [Result,varargout]=SPDE_run(param,time,space,BM,SPDE,methods,varargin)
%% Evaluates for declared methods and given parameters the following SPDE
% 
% $$ du(t,x)=\frac{a}{2} (\partial_{xx} u)(t,x) dt + \sigma (\partial_{x} u)(t,x) dW_t $$
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
mem_eff=0;
if ~isempty(varargin)
    for k=1:2:length(varargin)
        switch varargin{k}
            case 'memory'
                mem_eff=varargin{k+1};
        end
    end
end
mem_exact_enable=0;
mem_euler_enable=0;
mem_magnus_enable=1;
if mem_eff>=param.M
    mem_eff=param.M;
end
%
for k=1:1:length(methods)
    method=methods{k};
    fprintf('Calculating %s\n',method{1});
    switch method{1}
        case 'exact'
            if mem_eff>=1 && mem_exact_enable
                sstart=1;
                send=sstart+mem_eff-1;
                X=zeros(param.d,param.d,param.N,param.M);
                tic;
                while(send<=param.M)
                    w=sstart:1:send;
                    X(:,:,:,w)=...
                        SPDE_exact(...
                            SPDE.a,...
                            SPDE.sigma,...
                            time.t,...
                            space.x,...
                            BM.W(1,1,:,w),...
                            method{2}{:});
                    sstart=send+1;
                    send=sstart+mem_eff-1;
                end
                if sstart< param.M
                    w=sstart:1:param.M;
                    X(:,:,:,w)=...
                        SPDE_exact(...
                            SPDE.a,...
                            SPDE.sigma,...
                            time.t,...
                            space.x,...
                            BM.W(1,1,:,w),...
                            method{2}{:});
                end
                Result.exact.X=X; clear X;
                Result.exact.ctime.total=toc;
                fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            else
                tic;
                [Result.exact.X,Result.exact.failtimes]=...
                    SPDE_exact(...
                        SPDE.a,...
                        SPDE.sigma,...
                        time.t,...
                        space.x,...
                        BM.W,...
                        method{2}{:});
                Result.exact.ctime.total=toc;
                fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            end
        case 'euler'
            if mem_eff>=1 && mem_euler_enable
                sstart=1;
                send=sstart+mem_eff-1;
                X=zeros(param.d,param.d,param.N,param.M);
                tic;
                while(send<=param.M)
                    w=sstart:1:send;
                    X(:,:,:,w)=...
                        SPDE_euler(...
                            time,...
                            BM.W(1,1,:,w),...
                            BM.dW(1,1,:,w),...
                            SPDE.A,...
                            SPDE.B,...
                            method{2}{:});
                    sstart=send+1;
                    send=sstart+mem_eff-1;
                end
                if sstart< param.M
                    w=sstart:1:param.M;
                    X(:,:,:,w)=...
                        SPDE_euler(...
                            time,...
                            BM.W(1,1,:,w),...
                            BM.dW(1,1,:,w),...
                            SPDE.A,...
                            SPDE.B,...
                            method{2}{:});
                end
                Result.euler.X=X; clear X;
                Result.euler.ctime.total=toc;
                fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            else
                tic;
                Result.euler.X=...
                    SPDE_euler(...
                        time,...
                        BM.W,...
                        BM.dW,...
                        SPDE.A,...
                        SPDE.B,...
                        method{2}{:});
                Result.euler.ctime.total=toc;
                fprintf('Elapsed time %g\n',Result.(method{1}).ctime.total);
            end
        case {'m1','m2','m3'}
            m=method{1};
            fprintf('Calculating log\n');
%             tic;
%             Result.(method{1}).u=...
%                 SPDE_magnus_u0(time,BM.W,full(SPDE.A.sparse),full(SPDE.B.sparse),u0.u,str2num(m(end)),method{2}{:});
%             toc;
            if mem_eff>=1 && mem_magnus_enable
                sstart=1;
                send=sstart+mem_eff-1;
                X=zeros(param.d,param.d,param.N,param.M);
                logtime=0;
                exptime=0;
                while(send<=param.M)
%                 for w=1:1:param.M
                    w=sstart:1:send;
                    O=zeros(param.d,param.d,param.N,length(w));
                    fprintf('Trajectory %i up to %i from %i\n',w(1),w(end),param.M);
                    tic;
                    O(:,:,:,:)=...
                        SPDE_magnus(...
                        time.t,...
                        BM.W(:,:,:,w),...
                        full(SPDE.A.sparse),...
                        full(SPDE.B.sparse),...
                        str2num(m(end)),...
                        method{2}{:});
                    logtime=logtime+toc;
                    tic;
                    X(:,:,:,w)=m_exp(...
                    O,...
                    'type','matlab',...
                    'Comp Device','cpuReshape');
                    exptime=exptime+toc;
                    sstart=send+1;
                    send=sstart+mem_eff-1;
                end
                if sstart< param.M
                    w=sstart:1:param.M;
                    O=zeros(param.d,param.d,param.N,length(w));
                    fprintf('Trajectory %i up to %i from %i\n',w(1),w(end),param.M);
                    tic;
                    O(:,:,:,:)=...
                        SPDE_magnus(...
                        time.t,...
                        BM.W(:,:,:,w),...
                        full(SPDE.A.sparse),...
                        full(SPDE.B.sparse),...
                        str2num(m(end)),...
                        method{2}{:});
                    logtime=logtime+toc;
                    tic;
                    X(:,:,:,w)=m_exp(...
                    O,...
                    'type','matlab',...
                    'Comp Device','cpuReshape');
                    exptime=exptime+toc;
                end
                clear O;
                Result.(method{1}).X=X; clear X;
                Result.(method{1}).ctime.log=logtime;
                Result.(method{1}).ctime.expm=exptime;
                Result.(method{1}).ctime.total=...
                    Result.(method{1}).ctime.expm+...
                    Result.(method{1}).ctime.log;
            else
                tic;
                O=...
                    SPDE_magnus(...
                        time.t,...
                        BM.W,...
                        full(SPDE.A.sparse),...
                        full(SPDE.B.sparse),...
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
            end
        otherwise
            disp('Unknown method');
    end
end
end