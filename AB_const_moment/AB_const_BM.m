function [W,varargout]=AB_const_BM(time,M,varargin)
%% Creates M paths of the Browian motion for given time
%% Input: 
% # (struct) time: containing the fields t and dt, see AB_const_initialize
% # (integer) M_fine: containing number of paths
% # (cell) varargin: not needed
%% Output:
% # (1 x 1 x N x M array) W: containing M_fine paths of a Browian motion
% # (cell) varargout: dW: contatining the differences in time
%%
%
    W=zeros(1,1,size(time.t,3),M);
    W(1,1,2:end,:)=sqrt(time.dt).*randn(1,size(time.t,3)-1,M);
    if nargout > 1
        varargout{1}=W(1,1,2:end,:);
    end
    W=cumsum(W,3);
end