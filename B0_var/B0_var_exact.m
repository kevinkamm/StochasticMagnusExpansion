function X=B0_var_exact(f11,f12,f22,timegrid,W,varargin)
    comp_device='cpu';
    N=size(W,3);
    M=size(W,4);
    f11=reshape(f11,[N 1]);
    f12=reshape(f12,[N 1]);
    f22=reshape(f22,[N 1]);
    W=reshape(W,[N M]);
    dW=diff(W,1);
    timegrid=reshape(timegrid,[N 1]);
    dt=diff(timegrid);
    if ~isempty(varargin)
        for kk=1:1:length(varargin)
            switch varargin{kk}
                case 'Comp Device'
                    comp_device=varargin{kk+1};
            end
        end
    end
    if strcmp(comp_device,'cpu')
        Z=zeros(N,M);
    else
        Z=gpuArray.zeros(N,M);
        dW=gpuArray(dW);
        dt=gpuArray(dt);
    end
    x11=exp(...
        s_int(f11)-...
        l_int(f11.^2./2));%pointwise exp 
    x22=exp(...
        s_int(f22)-...
        l_int(f22.^2./2));%pointwise exp 
    quotient=x22./x11;
    x12=x11.*...
        (s_int(quotient.*f12)-...
        l_int(quotient.*f11.*f12));
    X=[reshape(x11,[1 1 N M]) reshape(x12,[1 1 N M]);
        reshape(Z,[1 1 N M]) reshape(x22,[1 1 N M])];
    X=gather(X);
    function I=l_int(X)
        I=Z;
        if size(X,2)==1
            I(2:end,:)=repmat(cumsum(X(1:end-1,:).*dt,1),[1 M]);
        else
            I(2:end,:)=cumsum(X(1:end-1,:).*dt,1);
        end
    end
    function I=s_int(X)
        I=Z;
        I(2:end,:)=cumsum(X(1:end-1,:).*dW,1);
    end
end