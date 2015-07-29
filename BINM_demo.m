clear
% Choose a data set to test BINM.
% data_set = 'Collins';
data_set = 'Friedel';


lambda = 1;
rho = 0.001;
max_iter = 20;
switch data_set
    case 'Collins'
        
        
        % Test BINM using the Collins dataset ('.\data\Collins_Scores.txt').
        % The strengths of direct interactions will be written into file  'Collins_result.txt' in current folder.
        W_dir = BINM_main('.\data\Collins_Scores.txt', 'Collins_result.txt', lambda, rho, max_iter);
        
        
    case 'Friedel'
        
        
        % Test BINM using the Friedel dataset ('.\data\Friedel_Scores.txt').
        % The strengths of direct interactions will be written into file  'Friedel_result.txt' in current folder.
        W_dir = BINM_main('.\data\Friedel_Scores.txt', 'Friedel_result.txt', lambda, rho, max_iter);
        
end