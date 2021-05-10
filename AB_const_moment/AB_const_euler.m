function X=AB_const_euler(dt,dW,A,B,varargin)
%% Function for evaluating the euler method
%% Input: 
% # see test_AB_const, AB_const_param_initialize, AB_const_coeff_initialize
% # (cell) varargin: Comp Device: either 'cpu' or 'gpu'
%% Output:
% # (d x N x M array) u
%%
%
    comp_device='cpu';
    disp = 0;
    if isempty(varargin)==0
        for k=1:1:length(varargin)
            switch varargin{k}
                case 'Comp Device'
                    comp_device=varargin{k+1};
                case 'disp'
                    disp=varargin{k+1};
            end
        end
    end
    switch comp_device
        case 'cpu'
            X=zeros(size(A,1),size(A,1),size(dW,3)+1,size(dW,4));
            X(:,:,1,:)=...
                repmat(...
                    reshape(...
                        eye(size(A,1),size(A,1)),...
                        [size(A,1),size(A,1),1,1]...
                    ),...
                    [1 1 1 size(dW,4)]...
                );
        case 'gpu'
            X=gpuArray.zeros(size(A,1),size(A,1),size(dW,3)+1,size(dW,4));
            X(:,:,1,:)=gpuArray(...
                        repmat(...
                            reshape(...
                                eye(size(A,1),size(A,1)),...
                                [size(A,1),size(A,1),1,1]...
                            ),...
                            [1 1 1 size(dW,4)]...
                        )...
                    );
    end
    for i=1:1:size(dW,3)
        if disp
            tic;
        end
        X(:,:,i+1,:)=euler_step(X(:,:,i,:),...
            mult(X(:,:,i,:),B,comp_device),...
            mult(X(:,:,i,:),A,comp_device),...
            dt,...
            dW(1,1,i,:));
%         ctime=toc;
        if disp
            fprintf('Time for evaluation %d seconds in step %d/%d\n',...
                ctime,i,size(dW,3));
        end
    end
    X=gather(X);
end
function X=euler_step(X,B,A,dt,dW)
    X=X+...
        B.*dt+...
        A.*dW;
end
function Z=mult(X,A,comp_device)
    switch comp_device
        case {'gpu'}
            Z=pagefun(@mtimes,A,X); 
        case {'cpu'}
            Z=pagemtimes(A,X);
        otherwise
            error('Unknown device')
    end
end