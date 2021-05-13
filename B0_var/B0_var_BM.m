function [W_fine,varargout]=B0_var_BM(time,M_fine,varargin)
%% Creates M paths of the Browian motion for given time
%% Input: 
% # (struct) time: containing the fields t_fine and dt_fine, see B0_var_initialize
% # (integer) M_fine: containing number of paths
% # (cell) varargin: not needed
%% Output:
% # (1 x 1 x N x M array) W_fine: containing M_fine paths of a Browian motion
% # (cell) varargout: dW_fine: contatining the differences in time
%%
%
    W_fine=zeros(1,1,size(time.t_fine,3),M_fine);
    W_fine(1,1,2:end,:)=sqrt(time.dt_fine).*randn(1,size(time.t_fine,3)-1,M_fine);
    if nargout > 1
        varargout{1}=W_fine(1,1,2:end,:);
    end
    W_fine=cumsum(W_fine,3);
end