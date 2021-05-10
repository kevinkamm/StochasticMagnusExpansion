function [time,BM]=AB_const_param_initialize(param)
% Initializes time and Brownian motion
% Input:
% # see test_AB_const
% Output:
% # (struct) time: fields: t and dt
% # (1 x 1 x N x 1 array) time.t: containing the time line
% # (1 x 1 x 1 x 1 double) time.dt: containing the step size
% # (struct) BM: fields: W and dW
% # (1 x 1 x N x M array) BM.W: containing the Browian motion
% # (1 x 1 x N-1 x M array) BM.dW: containing the differences in time of W
%

    time.t=...
        reshape(...
            linspace(param.t0,param.T,param.N),...
            [1 1 param.N 1]...
        );
    time.dt=(time.t(end)-time.t(1))/(param.N-1);
    [BM.W, BM.dW]=AB_const_BM(time,param.M);
end