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
    if isempty(varargin)==0
        for k=1:1:length(varargin)
            switch varargin{k}
                case 'Comp Device'
                    comp_device=varargin{k+1};
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
        X(:,:,i+1,:)=euler_step(X(:,:,i,:),...
            mult(X(:,:,i,:),B,comp_device),...
            mult(X(:,:,i,:),A,comp_device),...
            dt,...
            dW(1,1,i,:));
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
        case 'gpu'
            Z=pagefun(@mtimes,A,X); 
        otherwise
            Z=zeros(size(X));
            for j=1:1:size(X,4)
                Z(:,:,1,j)=A*X(:,:,1,j);
            end
    end
end
% function Z=mult(X,lower,main,upper,comp_device)
%     Z=toeplitz_tri_mmult(lower,main,upper,X,comp_device);
% end
% function T=toeplitz_tri_mmult(lower,main,upper,X,comp_device)
%     if strcmp(comp_device,'cpu')
%         Z=zeros(size(X));
%     else
%         Z=gpuArray.zeros(size(X));
%     end
%     %lower diag
%     temp=Z;
%     temp(2:end,1,:)=lower.*X(1:end-1,1,:);
%     temp(2:end,2:end,:)=lower.*X(1:end-1,2:end,:);
%     T=temp;
%     %main diag
%     T=T+main.*X;
%     %upper diag
%     temp=Z;
%     temp(1:end-1,1:end-1,:)=upper.*X(2:end,1:end-1,:);
%     temp(1:end-1,end,:)=upper.*X(2:end,end,:);
%     T=T+temp;
% end