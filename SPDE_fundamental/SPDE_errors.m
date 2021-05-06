function [Errors,varargout]=SPDE_errors(Result,time,errors,varargin)
%% Evaluates Errors of types: (absolute) L^p and (relative) L^1 for given Results
%% Input:
% * see test_AB_const and AB_const_run
%% Output:
% * (struct) Errors: Errors.(ip).(reference).(method).*
%%
% # ip: ('i%d_%d') where %d corresponds e.g. to the path i=151, j=51
% # reference: ('name') of reference method for error analysis
% # method: ('name') of method compared to reference for error analysis
% # * : either ('l%d_error') with default %d=2 or ('l1_error_rel')
% # (double) ('l%d_error') containing the L^p([0,T]\times \Omega,\bf(R)) error
% # (struct) ('l1_error_rel') containing fields: min, max, mean
%%
% internal paramters
p=2; % absolute error in L^p([0,T]\times\Omega;\bf(R)) 
cutlvl=0.1; % cut off between zero and cutlvl for relative errors
T=1; % set to terminal time to change L^p norm to a mean L^p norm
%%
%
if ~isempty(varargin)
    for k=1:1:length(varargin)
        switch varargin{k}
            case 'p'
                p=varargin{k+1};
            case 'cutlvl'
                cutlvl=varargin{k+1};
            case 'T'
                T=varargin{k+1};
        end
    end
end
varargout{1}=p;
varargout{2}=cutlvl;
varargout{3}=T;
lp_error=sprintf('l%i_error',p);
for k=1:1:length(errors)
    error=errors{k};
    pos1=error{1};
    pos2=error{2};
    fprintf('Displaying Errors for X(%i,%i,:,:)\n',pos1,pos2);
    reference=error{3};
    methods=error{4};
    X_ref=Result.(reference).X(pos1,pos2,:,:);
    fprintf('\t Reference method: %s\n',reference);
    for i=1:1:length(methods)
        fprintf('\t\t %s:\n',methods{i});
        X=Result.(methods{i}).X(pos1,pos2,:,:);
        ip=sprintf('p%i_%i',pos1,pos2);
        Errors.(ip).(reference).(methods{i}).(lp_error)=...
            lp_error_abs(X_ref,X,time,p)./T;
        fprintf('\t\t\t L_%i-error absolute:\n\t\t\t\t %g\n',...
            p,Errors.(ip).(reference).(methods{i}).(lp_error));
        [err_min,err_max,err_mean]=...
                    l1_error_rel(X_ref,X,time.dt,cutlvl);
        Errors.(ip).(reference).(methods{i}).l1_error_rel.min=...
            err_min;
        Errors.(ip).(reference).(methods{i}).l1_error_rel.max=...
            err_max;
        Errors.(ip).(reference).(methods{i}).l1_error_rel.mean=...
            err_mean./T;
        fprintf('\t\t\t L_1-error relative:\n');
        fprintf('\t\t\t\t min: %g\n',...
            Errors.(ip).(reference).(methods{i}).l1_error_rel.min);
        fprintf('\t\t\t\t max: %g\n',...
            Errors.(ip).(reference).(methods{i}).l1_error_rel.max);
        fprintf('\t\t\t\t mean: %g\n',...
            Errors.(ip).(reference).(methods{i}).l1_error_rel.mean);
    end
end
end
function err=lp_error_abs(X_ref,X,time,p)
    err=sum(mean(abs(X_ref-X).^p,4),3).*time.dt;
    err=err.^(1/p);
end
function [err_min,err_max,err_mean]=...
                        l1_error_rel(X_true,X_approx,dt,cutlvl)
    X_true=reshape(X_true,[size(X_true,3),size(X_true,4)]);
    X_approx=reshape(X_approx,[size(X_approx,3),size(X_approx,4)]);
    t_error_rel=-ones(1,size(X_true,1));
    for i=1:1:size(X_true,1)
        ind_cut=abs(X_true(i,:))>=cutlvl;
        if sum(ind_cut)>0
            t_error_rel(i)=...
                mean(abs(X_true(i,ind_cut)-X_approx(i,ind_cut))./...
                abs(X_true(i,ind_cut)),2);
        end
    end
    err_min=min(t_error_rel(t_error_rel>=0));
    if size(err_min,1)==0 || size(err_min,2)==0
        err_min=0;
    end
    err_max=max(t_error_rel);
    err_mean=sum(t_error_rel(t_error_rel>=0).*dt(1));
end