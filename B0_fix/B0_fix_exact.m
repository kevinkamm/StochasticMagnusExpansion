function X=B0_fix_exact(timegrid,W,varargin)
    N=size(W,3);
    M=size(W,4);
    W=reshape(W,[N M]);
    timegrid=reshape(timegrid,[N 1]);
    dt=diff(timegrid);
    dW=diff(W,1);
    comp_device='cpu';
    if ~isempty(varargin)
        for kk=1:1:length(varargin)
            switch varargin{kk}
                case 'Comp Device'
                    comp_device=varargin{kk+1};
            end
        end
    end
    if strcmp(comp_device,'cpu')
        Z=zeros(size(W));
    else
        Z=gpuArray.zeros(size(W));
        W=gpuArray(W);
        dW=gpuArray(dW);
        dt=gpuArray(dt);
        timegrid=gpuArray(timegrid);
    end   
    x11=exp(2.*(W-timegrid));
    x22=exp(-W-timegrid./2);
    integrand=...
        timegrid.*...
        exp((3/2).*timegrid-3.*W);
    x12=x11.*(s_int(integrand)-2.*l_int(integrand));
    X=[reshape(x11,[1 1 N M]) reshape(x12,[1 1 N M]);
        reshape(zeros(N,M),[1 1 N M]) reshape(x22,[1 1 N M])];
    X=gather(X);
    function I=l_int(X)
        I=Z;
        I(2:end,:)=cumsum(X(1:end-1,:).*repmat(dt,[1 M]),1);
    end
    function I=s_int(X)
        I=Z;
        I(2:end,:)=cumsum(X(1:end-1,:).*dW,1);
    end
end