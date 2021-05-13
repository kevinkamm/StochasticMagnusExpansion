function [time,BM,ind_t,ind_p]=AB_const_param_initialize(param)
%% Initializes time and Brownian motion
%% Input:
% # see test_AB_const
%% Output:
% # (struct) time: fields: t, t_fine and dt, dt_fine
% # (1 x 1 x N x 1 array) time.t: containing the time line
% # (1 x 1 x 1 x 1 double) time.dt: containing the step size
% # (1 x 1 x N_fine x 1 array) time.t_fine: containing the time line
% # (1 x 1 x 1 x 1 double) time.dt_fine: containing the step size
% # (struct) BM: fields: W, W_fine and dW, dW_fine
% # (1 x 1 x N x M array) BM.W: containing the Browian motion
% # (1 x 1 x N-1 x M array) BM.dW: containing the differences in time of W
% # (1 x 1 x N_fine x M_fine array) BM.W_fine: containing the Browian motion
% # (1 x 1 x N_fine-1 x M_fine array) BM.dW_fine: containing the differences 
% in time
% # (1 x N logical array) ind_t: containing conversation indizes from fine
% to normal for time axis, i.e. time.t=time.t_fine(ind_t)
% # (1 x M logical array) ind_p: containing conversation indizes from fine
% to normal for path axis, i.e. BM.W=BM.W_fine(ind_t,ind_m)
%% 
%
    time.t_fine=...
        reshape(...
            linspace(param.t0,param.T,param.N_fine),...
            [1 1 param.N_fine 1]...
        );
    time.dt_fine=(param.T-param.t0)/(param.N_fine-1);
%     ind_t=...
%     [1:1:param.N].*floor(param.N_fine/param.N)-...
%         (floor(param.N_fine/param.N)-1);
    ind_t=[1:1:param.N];
    ind_t(2:1:end)=ind_t(1:1:end-1).*floor((param.N_fine-1)/(param.N-1))+1;
    time.t=reshape(time.t_fine(:,:,ind_t,:),[1 1 param.N 1]);
    dt=diff(time.t);
    time.dt=dt(1); clear dt;
    [BM.W_fine, BM.dW_fine]=AB_const_BM(time,param.M_fine);
    ind_p=[1:1:param.M];
    BM.W=reshape(BM.W_fine(:,:,ind_t,ind_p),[1 1 param.N param.M]);
    BM.dW=diff(BM.W,1,3);
end