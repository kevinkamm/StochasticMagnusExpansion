function X=SPDE_euler(time,W,dW,A,B,varargin)
%% Function for evaluating the euler method
%%
% * Input: 
%%
% # see test_SPDE, SPDE_param_initialize, SPDE_coeff_initialize
% # (cell) varargin: Comp Device: either 'cpu' or 'gpu'
%%
% * Output:
%%
% # (d x d x N x M array) X
%%
%
    comp_device='cpu';
    d=size(A.sparse,1);
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
            X=zeros(d,d,size(W,3),size(W,4));
            X(:,:,1,:)=repmat(reshape(eye(d),[d d 1 1]),[1,1,1, size(W,4)]);
        case 'gpu'
            X=gpuArray.zeros(d,d,size(W,3),size(W,4));
            X(:,:,1,:)=gpuArray(repmat(reshape(eye(d),[d d 1 1]),[1,1,1, size(W,4)]));
    end
    dW=reshape(dW,[1 1 size(W,3)-1 size(W,4)]);
    for i=1:1:size(time.t,3)-1
        X(:,:,i+1,:)=euler_step(X(:,:,i,:),...
            mult(X(:,:,i,:),B.lower,B.main,B.upper,comp_device),...
            mult(X(:,:,i,:),A.lower,A.main,A.upper,comp_device),...
            time.dt,...
            dW(1,1,i,:));
    end
    X=gather(X);
end
function X=euler_step(X,B,A,dt,dW)
    X=X+...
        B.*dt+...
        A.*dW;
end
function Z=mult(X,lower,main,upper,comp_device)
    Z=toeplitz_tri_mmult(lower,main,upper,X,comp_device);
end
function T=toeplitz_tri_mmult(lower,main,upper,X,comp_device)
    if strcmp(comp_device,'cpu')
        Z=zeros(size(X));
    else
        Z=gpuArray.zeros(size(X));
    end
    %lower diag
    temp=Z;
    temp(2:end,1,1,:)=lower.*X(1:end-1,1,1,:);
    temp(2:end,2:end,1,:)=lower.*X(1:end-1,2:end,1,:);
    T=temp;
    %main diag
    T=T+main.*X;
    %upper diag
    temp=Z;
    temp(1:end-1,1:end-1,1,:)=upper.*X(2:end,1:end-1,1,:);
    temp(1:end-1,end,1,:)=upper.*X(2:end,end,1,:);
    T=T+temp;
end