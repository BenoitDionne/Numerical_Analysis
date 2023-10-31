% Usage:  arrayToLaTeXN(filename, matrix, cols, nbr, flag, caption, label,
%                     firsthead, secondhead)
%
% Produce LaTeX tables or longtables from an ( m x n )-matrix, where
% the columns used from the matrix are given in the vector cols.
% - 'string' : filename is the name of the file where the output will
%   be written.
% - The columns and rows of the LaTeX table will come from matrix.
% - [n1, ..., nk] : cols are the columns of the given matrix that will
%   be used to generate the LaTeX table.  The columns numbers given in
%   cols must be between 1 and n inclusively.  No column can be repeated.
% - number : The number of lines is split to form nbr tables side by
%   side to fill the width of a page.
% - If flag = 0, the table environment is used.
%   If flag <> 0, the longtable environment is used.
% - ['head1','head2',\ldots,'headk'] : The heading of each column for the
%   top of the table should be given firsthead.  the header for the
%   continuation of the table in a longtable environment should be given
%   in secondhead.  If secondhead is empty, the first head is used for the
%   continuation of the table.
% - 'string' : The caption is the caption for the LaTeX table.
% - 'string' : The label is the label for the LaTeX table.
%
function arrayToLaTeX(filename, matrix, cols, nbr, flag, caption, label, firsthead, secondhead)
    lth = size(matrix,1);
    Nrow = ceil(lth/nbr);
    Frow = floor(lth/nbr);
    Mcol = length(cols);
    Fcol = rem(lth,nbr);

    fid = fopen(filename,'w');

    if ( size(firsthead) ~= size(secondhead) )
        secondhead = firsthead;
    end

    if ( flag ~= 0 )
        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\setlongtables\n');
        fprintf(fid, '\\begin{longtable}{|');
    else
        fprintf(fid, '\\begin{table}\n');
        fprintf(fid, '\\caption{%s}\n', caption);
        fprintf(fid, '\\label{%s}\n', label);
        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{|');
    end

    for i=1:1:((Mcol+1)*nbr-1)
        fprintf(fid, 'c|');
    end
    fprintf(fid, '}\n');

    if ( flag ~= 0 )
        fprintf(fid, '\\caption{%s}\n', caption);
        fprintf(fid, '\\label{%s} \\\\ \n', label);
    end

    for i=1:1:nbr
        fprintf(fid, '\\cline{%d-%d}', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
    end
    fprintf(fid, '\n');

    for j=1:1:(nbr-1)
        for k=1:1:Mcol
            fprintf(fid, '%s & ', firsthead(k));
        end
        fprintf(fid, ' & ');
    end
    for k=1:1:(Mcol-1)
        fprintf(fid, '%s & ', firsthead(k));
    end
    fprintf(fid, '%s \\\\ \n', firsthead(Mcol));

    for i=1:1:nbr
        fprintf(fid, '\\cline{%d-%d}', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
    end
    fprintf(fid, '\n');

    if ( flag ~= 0 )
        fprintf(fid, '\\endfirsthead\n');

        fprintf(fid, '\\multicolumn{%d}{l}{\\small Continued from previous page.} \\\\ \n',(Mcol+1)*nbr-1);
        for i=1:1:nbr
            fprintf(fid, '\\cline{%d-%d} ', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
        end
        fprintf(fid, '\n');

        for j=1:1:(nbr-1)
            for k=1:1:Mcol
                fprintf(fid, '%s & ', secondhead(k));
            end
            fprintf(fid, ' & ');
        end
        for k=1:1:(Mcol-1)
            fprintf(fid, '%s & ', secondhead(k));
        end
        fprintf(fid, '%s \\\\ \n', secondhead(Mcol));
        for i=1:1:nbr
            fprintf(fid, '\\cline{%d-%d} ', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
        end
        fprintf(fid, '\n');
        fprintf(fid, '\\endhead\n');

        for i=1:1:nbr
            fprintf(fid, '\\cline{%d-%d}', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
        end
        fprintf(fid, '\n');
        fprintf(fid, '\\multicolumn{%d}{r}{\\small Continued on next page.} \\\\ \n',(Mcol+1)*nbr-1);
        fprintf(fid, '\\endfoot\n');

        for i=1:1:nbr
            fprintf(fid, '\\cline{%d-%d}', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
        end
        fprintf(fid, '\n');
        fprintf(fid, '\\endlastfoot\n');
    end

    for i=1:1:Nrow-1
        for j=1:1:nbr
            if j == 1
                N = i;
            elseif j>1 && j<(Fcol+2)
                N = N + Nrow;
            else
                N = N + Frow;
        end
        for k=1:1:(Mcol-1)
            fprintf(fid, '$%.8f$ & ', matrix(N,cols(k)));
	end
	fprintf(fid, '$%.8f$', matrix(N,cols(Mcol)));
	if ( j < nbr )
            fprintf(fid, ' & & ');
	end
     end
     fprintf(fid, ' \\\\ \n');
    end

    for j=1:1:nbr
        if ( Fcol<1 || j <= Fcol )
            N = j*Nrow;
            for k=1:1:(Mcol-1)
                fprintf(fid, '$%.8f$ & ', matrix(N,cols(k)));
            end
            fprintf(fid, '$%.8f$', matrix(N,cols(Mcol)));
        else
            for k=1:1:(Mcol-1)
                fprintf(fid, ' & ');
            end
        end
        if ( j < nbr )
            fprintf(fid, ' & & ');
        end
    end
    fprintf(fid, ' \\\\ \n');

    if ( flag ~= 0 )
        fprintf(fid, '\\end{longtable}\n');
        fprintf(fid, '\\end{center}\n');
    else
        for i=1:1:nbr
            fprintf(fid, '\\cline{%d-%d}', (i-1)*(Mcol+1)+1, i*(Mcol+1)-1);
        end
        fprintf(fid, '\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        fprintf(fid, '\\end{table}\n');
    end

    fclose(fid);

end

