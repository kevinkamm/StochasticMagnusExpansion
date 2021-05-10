function B0_fix_output(param,Result,output,varargin)
%% Creates output via LaTeX-template template.tex
% * available types: pdf
%% Input:
% * see test_AB_consts and AB_const_run
% * varargin: 'Errors', 'Plots', 'Surfaces', see AB_const_errors, 
% AB_const_plots, AB_const_surfaces
%%
%
    for kk=1:1:length(output)
        switch output{kk}
            case 'pdf'
                B0_fix_pdf(param,Result,varargin)
            otherwise
                disp('Output type unknown.');
        end
    end
end
function B0_fix_pdf(param,Result,varargin)
    pic_type='eps';
    save_param='epsc';
    root=[pwd,'\','Pdf'];
    fileName=...
            sprintf('B0_fix_T%1.3g_d%i_N%i_M%i',...
            param.T,param.d,param.N,param.M);
    tempPath=[root,'\','temp'];
    if exist(tempPath)==7
        rmdir(tempPath,'s');
    end
    mkdir(tempPath);
    inputPath=[root,'\','temp','\','input','.','tex'];
    delFile(inputPath);
    in_file=fopen(inputPath,'a');
    fprintf(in_file,...
        '\\section{Magnus expansion for $A_t\\in \\mathbb{R}^{2,2}$, $B=0$ deterministic}\n');
    filePath=pdf_Model();
    fprintf(in_file,...
        '\t\\input{%s}\n',changeSlash(filePath));
    fprintf(in_file,...
        '\\subsection{Parameters}\n');
    filePath=pdf_Paramters();
    fprintf(in_file,...
        '\t\\input{%s}\n',changeSlash(filePath));
    fprintf(in_file,...
        '\\subsection{Computational Times}\n');
    filePath=pdf_comp_times();
    fprintf(in_file,...
        '\t\\input{%s}\n',changeSlash(filePath));
    for kk=1:2:length(varargin{1})
        switch varargin{1}{kk}
            case 'Errors'
                fprintf(in_file,...
                    '\\subsection{Errors}\n');
                filePath=pdf_errors(varargin{1}{kk+1});
                pdf_errors_cmd(varargin{1}{kk+1});
                fprintf(in_file,...
                    '\t\\input{%s}\n',changeSlash(filePath));
            case 'Plots'
                fprintf(in_file,...
                    '\\subsection{Plots}\n');
                filePath=pdf_plots(varargin{1}{kk+1});
                fprintf(in_file,...
                    '\t\\input{%s}\n',changeSlash(filePath));
            case 'Surfaces'
                fprintf(in_file,...
                    '\\subsection{Surfaces}\n');
                filePath=pdf_surfaces(varargin{1}{kk+1});
                fprintf(in_file,...
                    '\t\\input{%s}\n',changeSlash(filePath));
            case 'Error_Plots'
                fprintf(in_file,...
                    '\\subsection{Error Plots}\n');
                filePath=pdf_error_plots(varargin{1}{kk+1});
                fprintf(in_file,...
                    '\t\\input{%s}\n',changeSlash(filePath));
            otherwise
                disp('Unknown argument.');
        end
    end
    fclose(in_file);
    templatePath=[pwd, '\' ,'template', '.','tex'];
    outputPath=[pwd,'\','template','.','pdf'];
%     tempPath=[root,'\','temp'];
    copyPath=[root,'\',fileName];
    copyPathFile=[root,'\',fileName,'\',fileName,'.','pdf'];
    if exist(copyPath)==7
        rmdir(copyPath,'s');
    end
    mkdir(copyPath);
%     delFile(copyPath);
    str1=sprintf('pdflatex %s',templatePath);
    system(str1);
    copyfile(tempPath,copyPath);
    copyfile(outputPath,copyPathFile);
    function filePath=pdf_Model()
        filePath=[root,'\','temp','\','model','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        fprintf(file,...
            'We will concern ourselves with the following SDE:\n');
        fprintf(file,...
            '\\begin{align*}\n');
        fprintf(file,...
            '\t dX_t = A_t X_t dW_t\n');
        fprintf(file,...
            '\\end{align*}\n');
        fprintf(file,...
            'with\n');
        str=...
            ['\\begin{align*}\n'...
            'A_t=\n'...
            '\\left[\n'...
            '\\begin{array}[c]{cc}\n'...
            'f^{11}(t) & f^{12}(t)\\\\\n'...
            '0 & f^{22}(t)\n'...
            '\\end{array}\n'...
            '\\right]\n'...
            '\\end{align*}\n'];
        fprintf(file,...
            str);
            fprintf(file,...
                'Here we chose:');
            str=...
            ['\\begin{align*}\n'...
            'f^{11}(t)\\equiv 2,\n'...
            'f^{12}(t)= t,\n'...
            'f^{22}(t)\\equiv -1,\n'...
            '\\end{align*}\n'];
            fprintf(file,...
                str);
        fclose(file);
    end
    function filePath=pdf_Paramters()
        filePath=[root,'\','temp','\','param','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        params=fieldnames(param);
        fprintf(file,...
            '\\begin{tabular}{@{}*{2}{c}@{}}\n');
        fprintf(file,...
            '\\text{\\textbf{Parameter}} & \\text{\\textbf{value}}\\\\\n');
        fprintf(file,...
            '\\toprule\\\\\n');
        for k=1:1:length(params)
            switch params{k}
                case 't0'
                    fprintf(file,...
                        '$t_0$ & $%g$\\\\\n',param.t0);
                case 'T'
                    fprintf(file,...
                        '$T$ & $%g$\\\\\n',param.T);
                case 'N'
                    fprintf(file,...
                        '$N$ & $%i$\\\\\n',param.N);
                case 'M'
                    fprintf(file,...
                        '$M$ & $%i$\\\\\n',param.M);
                case 'd'
                    fprintf(file,...
                        '$d$ & $%i$\\\\\n',param.d);
                case 'N_fine'
                    fprintf(file,...
                        '\\verb+N_fine+ & $%g$\\\\\n',param.N_fine);
                case 'M_fine'
                    fprintf(file,...
                        '\\verb+M_fine+ & $%g$\\\\\n',param.M_fine);
            end
        end
        fprintf(file,...
            '\\end{tabular}\n');
        fclose(file);
    end
    function filePath=pdf_comp_times()
        filePath=[root,'\','temp','\','ctimes','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        methods=fieldnames(Result);
        fprintf(file,...
            '\\begin{tabular}{@{}*{4}{c}@{}}\n');
        fprintf(file,...
            '\\text{\\textbf{Method}} &');
        fprintf(file,...
            '\\text{\\textbf{Log}} &');
        fprintf(file,...
            '\\text{\\textbf{Matrix Exp}} &');
        fprintf(file,...
            '\\text{\\textbf{Total}}\\\\\n');
        fprintf(file,...
            '\\toprule\\\\\n');
        for k=1:1:length(methods)
            fprintf(file,...
                '\\text{%s} & ',...
                methods{k});
            if isfield(Result.(methods{k}).ctime,'log')
                fprintf(file,...
                '%g & ',...
                Result.(methods{k}).ctime.log);
            else
                fprintf(file,...
                    '0 & ');
            end
            if isfield(Result.(methods{k}).ctime,'expm')
                fprintf(file,...
                '%g & ',...
                Result.(methods{k}).ctime.expm);
            else
                fprintf(file,...
                    '0 & ');
            end
            fprintf(file,...
                '%g \\\\\n',...
                Result.(methods{k}).ctime.total);
        end
        fprintf(file,...
            '\\end{tabular}\n');
        fclose(file);
    end
    function filePath=pdf_errors(Errors)
        filePath=[root,'\','temp','\','errors','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        ijs=fieldnames(Errors);
        fprintf(file,...
                '\\begin{compactenum}\n');
        for n=1:1:length(ijs)
            ip=ijs{n};
            if ~strcmp(ip(1),'p')
                fprintf(file,...
                    '\\item Total Errors:\n');
                references=fieldnames(Errors.total);
                fprintf(file,...
                '\\begin{compactenum}\n');
                for i=1:1:length(references)
                    fprintf(file,...
                        '\\item Reference method: %s\\\\\n',references{i});
                    methods=fieldnames(Errors.total.(references{i}));
                    fprintf(file,...
                        '\\begin{tabular}{@{}*{2}{c}@{}}\n');
                    fprintf(file,...
                        '\\text{\\textbf{Method}} & \\text{$\\mathbb{E}[Err_{%2.3g}]$}\\\\\n\\toprule\n',param.T);
                    for j=1:1:length(methods)
                        fprintf(file,...
                            '%s &',methods{j});
                        fprintf(file,...
                            '$%3.3g\\,\\%%$ \\\\\n',Errors.total.(references{i}).(methods{j}).expectation*100);
                    end
                    fprintf(file,...
                    '\\end{tabular}\n');
                end
                fprintf(file,...
                '\\end{compactenum}\n');
            else
                ij=split(ip(2:end),'_');
                ii=ij{1};
                jj=ij{2};
                fprintf(file,...
                        '\\item Errors for $X(%s,%s,:,:)$:\n',ii,jj);
                references=fieldnames(Errors.(ip));
                fprintf(file,...
                    '\\begin{compactenum}\n');
                for i=1:1:length(references)
                    fprintf(file,...
                        '\\item Reference method: %s\\\\\n',references{i});
                    methods=fieldnames(Errors.(ip).(references{i}));
                    fprintf(file,...
                        '\\begin{tabular}{@{}*{%i}{c}@{}}\n',length(methods)+1);
                    fprintf(file,...
                        '\\text{\\textbf{Error}} &');
                    err=zeros(4,length(methods));
                    types=fieldnames(Errors.(ip).(references{i}).(methods{1}));
                    for j=1:1:length(methods)
                        if j<length(methods)
                            fprintf(file,...
                                '\\text{\\textbf{%s}} &',methods{j});
                        else
                            fprintf(file,...
                                '\\text{\\textbf{%s}} \\\\\n',methods{j});
                        end
                        err(1,j)=...
                            Errors.(ip).(references{i}).(methods{j}).(types{1});
                        err(2,j)=...
                            Errors.(ip).(references{i}).(methods{j}).(types{2}).min;
                        err(3,j)=...
                            Errors.(ip).(references{i}).(methods{j}).(types{2}).max;
                        err(4,j)=...
                            Errors.(ip).(references{i}).(methods{j}).(types{2}).mean;
                    end
                    fprintf(file,...
                        '\\toprule\\\\\n');
                    for i1=1:1:size(err,1)
                        switch i1
                            case 1
                                fprintf(file,...
                                    '(abs error) L%s &',types{1}(2));
                            case 2
                                fprintf(file,...
                                    '(rel error) min &');
                            case 3
                                fprintf(file,...
                                    '(rel error) max &');
                            case 4
                                fprintf(file,...
                                    '(rel error) mean &');
                            otherwise
                                fprintf(file,...
                                    ' &');
                        end
                        for i2=1:1:size(err,2)
                            fprintf(file,...
                                '$%g$ ',err(i1,i2));
                            if(i2==size(err,2))
                                fprintf(file,...
                                    '\\\\\n');
                            else
                                fprintf(file,...
                                    '&');
                            end
                        end
                    end
                    fprintf(file,...
                        '\\end{tabular}\n');
                end
                fprintf(file,...
                    '\\end{compactenum}\n');
            end
        end
        fprintf(file,...
                '\\end{compactenum}\n');
        fclose(file);
    end
    function filePath=pdf_errors_cmd(Errors)
        Name='total_err_tab';
        filePath=[root,'\','temp','\',Name,'.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        ijs=fieldnames(Errors);
        for n=1:1:length(ijs)
            ip=ijs{n};
            if ~strcmp(ip(1),'p')
                references=fieldnames(Errors.total);
                for i=1:1:length(references)
                    methods=fieldnames(Errors.total.(references{i}));
                    fprintf(file,...
                        '\\begin{tabular}{@{}*{1}{c}@{}}\n');
                    fprintf(file,...
                        '\\text{\\ $\\mathbb{E}[Err_{%2.3g}]$}\\\\\n\\toprule\n',param.T);
                    for j=1:1:length(methods)
                        fprintf(file,...
                            '\\ $%3.3g\\,\\%%$ \\\\\n',Errors.total.(references{i}).(methods{j}).expectation*100);
                    end
                    fprintf(file,...
                    '\\end{tabular}%%\n');
                end
            end
        end
        fclose(file);
    end
    function filePath=pdf_plots(Plots)
        filePath=[root,'\','temp','\','plots','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        for i=1:1:length(Plots)
            fig_name=...
                sprintf('plot_%i',i);
            pngPath = [root,'\','temp','\',fig_name,'.',pic_type];
            delFile(pngPath);
            saveas(Plots(i),pngPath,save_param);
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(pngPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end
    function filePath=pdf_error_plots(Error_Plots)
        filePath=[root,'\','temp','\','error_plots','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        for i=1:1:length(Error_Plots)
            fig_name=...
                sprintf('error_plot_%i',i);
            pngPath = [root,'\','temp','\',fig_name,'.',pic_type];
            delFile(pngPath);
            saveas(Error_Plots(i),pngPath,save_param);
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(pngPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end
    function filePath=pdf_surfaces(Surfaces)
        filePath=[root,'\','temp','\','surfaces','.','tex'];
        delFile(filePath);
        file = fopen(filePath,'a');
        for i=1:1:length(Surfaces)
            fig_name=...
                sprintf('surface_%i',i);
            pngPath = [root,'\','temp','\',fig_name,'.',pic_type];
            delFile(pngPath);
            saveas(Surfaces(i),pngPath,save_param);
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(pngPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end
end
function str=changeSlash(str)
    for i=1:1:length(str)
        if strcmp(str(i),'\')
            str(i)='/';
        end
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function matrix2Latex(file,M)
    [n,m]=size(M);
    fprintf(file,...
        '\\left[\\begin{array}[c]{*{%i}{c}}\n',n);
    for i=1:1:n
        for j=1:1:m
            fprintf(file,...
                '%g ',M(i,j));
            if(j==m)
                fprintf(file,...
                    '\\\\\n');
            else
                fprintf(file,...
                    '&');
            end
        end
    end
    fprintf(file,...
        '\\end{array}\\right]\n');
end