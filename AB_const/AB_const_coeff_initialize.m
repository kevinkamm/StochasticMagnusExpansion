function [A,B]=AB_const_coeff_initialize(param,varargin)
%% Function for initializing the A and B
%% Input: 
% # see test_AB_const, AB_const_param_initialize
% # (cell) varargin: 'Example', 'load'
%% Output:
% # (d x d x 1 x 1 array) A: containing coefficient matrix A
% # (d x d x 1 x 1 array) B: containing coefficient matrix B
%%
%
    Example='default';
    reload=1;
    if ~isempty(varargin)
        for kk=1:1:length(varargin)
            switch varargin{kk}
                case 'Example'
                    Example=varargin{kk+1};
                case 'load'
                    reload=varargin{kk+1};
            end
        end
    end
    strA=...
        ['AB','\'...
        Example,'\',...
        sprintf('A_%i.mat',param.d)];
    strB=...
        ['AB','\',...
        Example,'\',...
        sprintf('B_%i.mat',param.d)];
    if reload==1
        [A,B]=loadExample();
    else
        [A,B]=createExample();
    end
    function [A,B]=loadExample()
        if exist(strA)~=2 || exist(strB)~=2
            [A,B]=createExample();
        else
            AA=load(strA);
            BB=load(strB);
            A=AA.A;
            B=BB.B;
            clear AA BB;
        end
    end
    function [A,B]=createExample()
        switch Example
            case 'default'
                A=randn(param.d);
                A=A./norm(A,2);
                B=randn(param.d);
                B=B./norm(B,2);
                save(strA,'A');
                save(strB,'B');
        end
    end
end