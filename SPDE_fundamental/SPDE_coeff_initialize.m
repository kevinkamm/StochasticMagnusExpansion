function SPDE=SPDE_coeff_initialize(param,space,SPDE)
%% Function for initializing the SPDE instance
%%
% * Input: 
%%
% # see test_SPDE, SPDE_param_initialize
%%
% * Output:
%%
% # (struct) SPDE: SPDE.([A,B]).([sparse,lower,main,upper])
%%
%
    SPDE.A.sparse=...
        spdiags(...
        [-ones(param.d,1) ones(param.d)],...
        [-1 0],...
        param.d,param.d).*(SPDE.sigma/space.dx);
    SPDE.B.sparse=...
        spdiags(...
        [ones(param.d,1)./2 -ones(param.d,1) ones(param.d)./2],...
        [-1 0 1],...
        param.d,param.d).*(SPDE.a/(space.dx)^2);
    SPDE.A.lower=-(SPDE.sigma/space.dx);
    SPDE.A.main=-SPDE.A.lower;
    SPDE.A.upper=0;
    SPDE.B.lower=SPDE.a/(2.*(space.dx)^2);
    SPDE.B.main=-(SPDE.a/(space.dx)^2);
    SPDE.B.upper=SPDE.B.lower;
end