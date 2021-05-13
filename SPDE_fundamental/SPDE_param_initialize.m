function [time,space,BM]=SPDE_param_initialize(param)
%% Initializes time, space, Brownian motion and initial datum
%% Input:
% # see test_SPDE
%% Output:
% # (struct) time: fields: t and dt
% # (1 x N x 1 array) time.t: containing the time line
% # (1 x 1 x 1 double) time.dt: containing the step size
% # (struct) space: fields: x and dx
% # (d+2 x 1 x 1 array) space.x: containg the spatial line
% # (1 x 1 x 1 double) space.dx: containg the step size
% # (1 x N x M array) BM: containing the Browian motion
% # (d x 1 x 1 array) u0: containing the initial datum
%% 
%
    time.t=reshape(linspace(param.t0,param.T,param.N),[1 1 param.N 1]);
    time.dt=(param.T-param.t0)/(param.N-1);
    space.x=reshape(linspace(param.xa,param.xb,param.d+2),[param.d+2 1 1]);
    space.dx=(param.xb-param.xa)/(param.d+1);
    [BM.W, BM.dW]=SPDE_BM(time,param.M);
end