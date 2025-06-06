function tab_print(matrix, filename)

% tab_print(matrix, filename)
%
% Prints "matrix" to the file "filename" in tab-delimited format.
%
% Tab-delimited data can be opened nicely in spreadsheet programs like
% Excel and SPSS.
%
% "matrix" should be any numerical Matlab matrix.
%
% "filename" should be a string, for instance "myvariable.txt"
%
% "filename" will be printed to the current Matlab directory.
% Alternatively, you can specify the print directory by specifying the full
% path in "filename," for instance "C:\MyFolder\myvariable.txt"
%
% If "filename" is not specified, "matrix" is printed in tab-delimited format
% to the Matlab prompt. This text can then be copy-pasted into a
% spreadsheet.

if ~exist('filename','var') || isempty(filename)
    fid = 1;
    saved_to_file = 0;
else
    fid = fopen(filename,'w');
    saved_to_file = 1;
end

fprintf(1,'\n');
for i = 1:size(matrix,1)
    fprintf(fid,'%f\t',matrix(i,:));
    if i < size(matrix,1)
        fprintf(fid,'\n');
    end
end

if saved_to_file
    fprintf(1,'great success');
    fclose(fid);    
end

fprintf(1,'\n\n\n');