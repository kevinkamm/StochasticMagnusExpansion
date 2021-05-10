function [A,varargout]=B0_var_coeff_initialize(param,time,varargin)
%% Function for initializing the A
%% Input: 
% # see test_B0_var, B0_var_param_initialize
% # (cell) varargin: 'Example'
%% Output:
% # (struct) A: containing the fields: f11,f12,f22,df11,df12,df22
%%
%
    A.f11=zeros(param.N_fine,1);
    A.f12=zeros(param.N_fine,1);
    A.f22=zeros(param.N_fine,1);
    A.df11=zeros(param.N_fine,1);
    A.df12=zeros(param.N_fine,1);
    A.df22=zeros(param.N_fine,1);
    t=reshape(time.t_fine,[param.N_fine,1]);
    Example='fix';
    if ~isempty(varargin)
        for kk=1:1:length(varargin)
            switch varargin{kk}
                case 'Example'
                    Example=varargin{kk+1};
            end
        end
    end
    switch Example
        case 'fix'
            fix();
            varargout{1}='fix';
        case 'scaled'
            fix();
            scaled();
            varargout{1}='scaled';
    end
    function fix()
        A.f11=2.*ones(param.N_fine,1);
        A.f12=t;
        A.f22=-ones(param.N_fine,1);
        A.df12=ones(param.N_fine,1);
    end
    function scaled()
        function sigma_t=spectral_radius(t)
%             sigma_t=zeros(param.N_fine,1);
            sigma_t=...
                sqrt(...
                    (t.^2+5)./2+...
                    sqrt(...
                        t.^4+...
                        10.*t.^2+...
                        9)./2);
        end
        sigma_t=feval(@(t) spectral_radius(t),t);
        A.f11=A.f11./sigma_t;
        A.f12=A.f12./sigma_t;
        A.f22=A.f22./sigma_t;
        function df=derivative1(t)
            numerator=...
                ((4.*t.^3+20.*t)./...
                (4.*sqrt(t.^4+10.*t.^2+9)))+...
                t;
            denominator=...
                2.*(...
                    (t.^2+5)./2+...
                    sqrt(...
                        t.^4+...
                        10.*t.^2+...
                        9)./2).^(3/2);
            df=(-numerator./denominator);
        end
        function df2=derivative2(t)
            denominator=...
                sqrt(...
                    (t.^2+5)./2+...
                    sqrt(...
                        t.^4+...
                        10.*t.^2+...
                        9)./2);
            df2=(ones(size(denominator))./denominator);
        end
        df=feval(@(t)derivative1(t),t);
        df2=feval(@(t)derivative2(t),t)+...
            t.*df;
        A.df11=2.*df;
        A.df12=df2;
        A.df22=-df;
    end
end