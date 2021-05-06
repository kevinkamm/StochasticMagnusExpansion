function [X,varargout]=SPDE_exact(a,sigma,timegrid,x,W,varargin)
%% Calculates exact solution.
% Please note, that this method does not work for small times and will 
% give results as NaN and Inf because e^x is beyond double precision for
% necessaray large inputs
    comp_device='cpu';
    m=size(x,1);
    N=size(timegrid,3);
    M=size(W,4);
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
            X=repmat(reshape(eye(m-2),[m-2 m-2 1 1]),[1 1 N M]);
        case 'gpu'
            X=gpuArray(repmat(reshape(eye(m-2),[m-2 m-2 1 1]),[1 1 N M]));
    end
%     maxcut=1;
    failtimes=logical(zeros(m-2,m-2,N,M));
%     xi=reshape(x(2:end-1),[1 m-2 1 1]);clear spacegrid;
    xi=reshape((x(1:1:end-1)+x(2:1:end))./2,[1 m-1 1 1]);
    t=reshape(timegrid(2:end),[1 1 N-1 1]);clear timegrid;
    xsW=reshape(x(2:end-1),[m-2 1 1 1])+reshape(sigma.*W(1,1,2:end,:),[1 1 N-1 M]);clear W;
    Ip=scaled_erf(xi);
    X(:,:,2:end,:)=...
        (Ip(:,1:1:end-1,:,:)-Ip(:,2:1:end,:,:))./2;
    indInf=isinf(X);
    indNaN=isnan(X);
    failtimes(indInf)=1;
    failtimes(indNaN)=1;
    varargout{1}=failtimes;
    X=gather(X);
    function Ip=scaled_erf(xi)
        Ip=erf((-xi+xsW)./...
            sqrt(2.*(a-sigma^2).*t));
%         ind=isinf(Ip);
%         failtimes(ind)=1;
%         ind=isnan(Ip);
%         failtimes(ind)=1;
    end
end