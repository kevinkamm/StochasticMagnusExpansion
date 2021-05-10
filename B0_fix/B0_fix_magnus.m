function O=B0_fix_magnus(timegrid,W,order,varargin)
% Requires B=0, A=[2 t; 0 -1] 
% Input:
%   timegrid
%   W
%   order
% Output:
%   order=1 then O d\times d \times N \times M tensor
%   order=2 then O d\times d \times N \times M tensor
%   O [d\times d \times N \times M,d\times d \times N \times M]
%       tensor
    N=size(W,3);
    M=size(W,4);
    W=reshape(W,[N M]);
    d=2;
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
        O=zeros(d,d,N,M);
    else
        O=gpuArray.zeros(d,d,N,M);
        W=gpuArray(W);
        dt=gpuArray(dt);
        timegrid=gpuArray(timegrid);
    end  
    function O1=firstorder()
        IW=l_int(W);
        O1=O;
        O1(1,1,:,:)=...
            2.*W;
        O1(1,2,:,:)=...
            timegrid.*W-IW;
        O1(2,2,:,:)=...
           -W;
        O1=gather(O1);
    end
    function O2=secondorder()
        IW=l_int(W);
        IW2=l_int(W.^2);
        O2=O;
        O2(1,1,:,:)=...
            repmat(-2.*timegrid,[1 M]);
        O2(2,2,:,:)=...
            repmat(timegrid./(-2),[1 M]);
        O2(1,2,:,:)=...
            timegrid.^2./(-4)+...
            (-3./2).*(W.*IW-IW2);
        O2=gather(O2);
    end
    function O3=thirdorder()
        IW=l_int(W);
        IW2=l_int(W.^2);
        IW3=l_int(W.^3);
%         IsW2=l_int(timegrid.*W.^2);
        IsW=l_int(timegrid.*W);
        O3=O;
%         O3(1,2,:,:)=...
%             (15./16).*IW.^2-...
%             (7./8).*(...
%                 timegrid.*IW2-...
%                 IsW2...
%             )-...
%             (1./24).*timegrid.^3+...
        O3(1,2,:,:)=...
            -(3./4).*W.^2.*IW+...
            (3./4).*timegrid.*IW-...
            (3./2).*IsW+...
            (9./4).*W.*IW2-...
            (3./2).*IW3+...
            (3./8).*timegrid.^2.*W;
        O3=gather(O3);
    end
    function I=l_int(W)
        I=reshape(O(1,1,:,:),[N M]);
        I(2:1:end,:)=cumsum(W(1:1:end-1,:).*dt,1); 
    end
    if order==1
        O=firstorder();
    elseif order==2
        O=firstorder()+secondorder();
    elseif order==3
        O=firstorder()+secondorder()+thirdorder();
    else
        O=[firstorder(),secondorder(),thirdorder()];
    end
end