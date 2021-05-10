function O=SPDE_magnus(timegrid,W,A,B,order,varargin)
%% Calculates the Logarithm O of the Magnus expansion
% Requires A,B to be deterministic and quadratic
%% Input:
% # (1 x 1 x N x 1 array) timegrid: containing the time line
% # (1 x 1 x N x M array) W: containing the Browian motion
% # (d x d x 1 x 1 array) A: containing the coefficient matrix A
% # (d x d x 1 x 1 array) B: containing the coefficient matrix B
% # (integer) order: containing the order of the Magnus expansion: 1,2,3
% # (cell) varargin: Comp Device: either 'cpu' or 'gpu'
%
%% Output:
% # (d x d x N x M array) O: containing the Magnus logarithm of given order
%%
%
d=size(A,1);
N=size(W,3);
M=size(W,4);
dt=reshape(diff(timegrid),[1 1 N-1 1]);
A=reshape(A,[d d 1 1]);
B=reshape(B,[d d 1 1]);
W=reshape(W,[1 1 N M]);
timegrid=reshape(timegrid,[1 1 N 1]);
comp_device='cpu';
if ~isempty(varargin)
    for k=1:1:length(varargin)
        switch varargin{k}
            case 'Comp Device'
                comp_device=varargin{k+1};
        end
    end
end
if strcmp(comp_device,'cpu')
    O=zeros(d,d,N,M);
    slicemtimes=str2func('bsxfun');
else
    O=gpuArray(zeros(d,d,N,M));
    A=gpuArray(A);
    B=gpuArray(B);
    dt=gpuArray(dt);
    W=gpuArray(W);
    timegrid=gpuArray(timegrid);
    slicemtimes=str2func('pagefun');
end
    function O1=firstorder()
%         O1=zeros(d,d,N,M);
        O1=slicemtimes(@mtimes,B,timegrid)+...
            A.*W;
        O1=gather(O1);
    end
    function O2=secondorder()
%         O2=zeros(d,d,N,M);
        BA=comm(B,A);
        I1=l_int(W);
        O2=...
            slicemtimes(@mtimes,-A^2./2,timegrid)+...
            BA.*I1-slicemtimes(@mtimes,BA,timegrid).*W./2;
        O2=gather(O2);
    end
    function O3=thirdorder()
        BA=comm(B,A);
        BAA=comm(BA,A);
        BAB=comm(BA,B);
        IsW=l_int(W.*timegrid);
        IW=l_int(W);
        IW2=l_int(W.^2);
        O3=...
           -slicemtimes(@mtimes,BAB./12,timegrid.^2).*W+...
          slicemtimes(@mtimes,BAA./12,timegrid).*...
            W.^2+...
          BAB.*...
            IsW-...
          BAA./2.*...
            W.*IW-...
          slicemtimes(@mtimes,BAB./2,timegrid).*...
            IW+...
          BAA./2.*...
            IW2;
        O3=gather(O3);
    end
    function I=l_int(W)
        I=O(1,1,:,:);
        I(1,1,2:1:N,:)=cumsum(W(1,1,1:1:end-1,:).*dt,3); 
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
function Z=comm(A,B)
    Z=A*B-B*A;
end