function E=m_exp(X,varargin)
%% Function for evaluating a matrix exponential e^X
%% Input: 
% # (d x d x N x M array) X containing the matrix for the matrix exp
% # (cell) varargin: type, Comp Device
% # possible types: matlab
% # possible Comp Devices: cpu
%% Output:
% # (d x d x N x M array) E=e^X
%%
%
    type='matlab';
    comp_device='cpu';
    if ~isempty(varargin)
        for k=1:1:length(varargin)
            switch varargin{k}
                case 'type'
                    type=varargin{k+1};
                case 'Comp Device'
                    comp_device=varargin{k+1};
            end
        end
    end
    switch type
        case 'matlab'
            switch comp_device
                case 'cpuReshape'
                    N=size(X,3);
                    M=size(X,4);
                    X=reshape(X,[size(X,1),size(X,2),N*M]);
                    E=zeros(size(X,1),size(X,2),size(X,3));
                    parfor i=1:size(X,3)
                        E(:,:,i)=expm(X(:,:,i));
                    end
                    E=reshape(E,[size(E,1),size(E,2),N,M]);
                case 'cpu'
                    E=zeros(size(X,1),size(X,2),size(X,3),size(X,4));
                    for j=1:1:size(X,4)
                        for i=1:1:size(X,3)
                            E(:,:,i,j)=expm(X(:,:,i,j));
                        end
                    end
                case 'gpu'
                    E=zeros(size(X,1),size(X,2),size(X,3),size(X,4));
                    for j=1:1:size(X,4)
                        for i=1:1:size(X,3)
                            E(:,:,i,j)=expm(X(:,:,i,j));
                        end
                    end
            end
    end
end