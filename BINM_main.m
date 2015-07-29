function W_dir = BINM_main(PPI_profie, output_file_name, lambda, rho, max_iter)
% GMFTP_main is the main function of our model. It reads data (AP-MS PPI network) from the text
% file (function: Data_Read.m) and identifies binary interactions using the model
% presented in the paper (function: BINM.m). It also writes the confidence scores of direct interactions into file 'output_file_name', where each line
% presents an interaction between two proteis with the weight of direct interaction. Example:
% YHR172W	YNL126W	0.445809

% Inputs:
%   PPI_profie: the input file name of the AP-MS PPI network, where each line
%   contains two proteins defining an interaction with the confidence score derived from AP-MS data.
%   Example: YAL001C    YBR123C	0.983749

%   output_file_name: the file into which GMFTP writes the confidence scores of direct interactions. Each line presents an interaction with strength of direct interaction.
%   Example:  YHR172W	YNL126W	0.445809
%   rho: the tolerance threshold of the stop criterion. The default value is 1e-5.
%	max_iter: the number of iterations limited in BINM. The default value is 400.

% Outputs:
%   W_dir: the estimator of direct interaction matrix.

if nargin < 5
    max_iter = 20;
end

if nargin < 4
    rho = 0.001;
end


if nargin < 3
    lambda = 1;
end


if nargin < 2
    output_file_name = 'result.txt';
end



% Read data from the text file and write it into matlab format
Data_set = Data_Read(PPI_profie);

% Identifying binary interactions using the propsed model.
[W_dir, ~] = BINM(Data_set.PPI, lambda, rho, max_iter);

% Write the confidence scores of direct interactions into file 'output_file_name'.
Result_Print(output_file_name, Data_set.Protein, W_dir);




function Data_set = Data_Read(PPI_profie)
fprintf('Reading data...')
fprintf('\n')
fid_ppi=fopen(PPI_profie);
temp_PPI=textscan(fid_ppi,'%s%s%f%*[^\n]','delimiter',{'\t',' '},'Headerlines',0);
fclose(fid_ppi);

Data_set.Protein = union(temp_PPI{1},temp_PPI{2});
[~,Locb_1] = ismember(temp_PPI{1}, Data_set.Protein);
[~,Locb_2] = ismember(temp_PPI{2}, Data_set.Protein);
Data_set.PPI = sparse(Locb_1,Locb_2,temp_PPI{3},length(Data_set.Protein),length(Data_set.Protein));
Data_set.PPI = Data_set.PPI + Data_set.PPI';
Data_set.PPI = Data_set.PPI - diag(diag(Data_set.PPI));


function   Result_Print(output_file_name, protein_list, S)

[I1,I2,C] = find(triu(S));
[Y, index] = sort(C, 'descend');
I1 = I1(index);
I2 = I2(index);
C = C(index);
fid = fopen(output_file_name,'w');
for k = 1:size(I1)
    fprintf(fid, '%s\t', cell2mat( protein_list(I1(k))) );
    fprintf(fid, '%s\t', cell2mat( protein_list(I2(k))) );
    fprintf(fid, '%f', C(k) );
    fprintf(fid, '\n');
end

fclose(fid);
fprintf(['The results have been written into file ', output_file_name])
fprintf('\n')
