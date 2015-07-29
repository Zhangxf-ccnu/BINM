function [W_dir, score] = BINM(W_obs, lambda, rho, max_iter)
% The core function of this paper. It estimates parameter W_{dir} from a given AP-MS PPI network with weighted adjacent matrix W_{obs}.

% Input:
%   W_obs: weighted adjacent matrix of AP-MS PPI network with size N by N.
%   rho: the tolerance threshold of the stop criterion. The default value is 1e-5.
%	max_iter: the number of iterations limited in BINM. The default value is 400.


% Outputs:
%   W_dir: the estimator of direct interaction matrix.
%   score: the value of objective function (3)

if nargin < 4
    max_iter = 20;
end

if nargin < 3
    rho = 0.001;
end


if nargin < 2
    lambda = 1;
end

if nargin < 1
    error('You need input W_obs');
end

W_obs = W_obs - diag(diag(W_obs)); % Set the diagonal terms to zero
W_dir = W_obs; %Initialization

loss_value_old = inf;

for i = 1:max_iter
i
    % Compute  W_obs_hat
    W_obs_hat = W_dir + (W_dir*W_dir);
    W_obs_hat = W_obs_hat - diag(diag(W_obs_hat));
    % Update W_dir according to Equation (4)
    W_dir = W_dir.* ( ( W_obs*W_dir' + W_dir'*W_obs+ W_obs) ./ (W_obs_hat*W_dir' + W_dir'*W_obs_hat+ W_obs_hat+ lambda*W_dir + eps) ).^(1/4);
    % Calculate the value of objective function (3)
    loss_value_new = function_eval(W_obs, W_dir, lambda);
     abs(loss_value_new - loss_value_old) / abs(loss_value_new) 
    if abs(loss_value_new - loss_value_old) / abs(loss_value_new) < rho
        break;
    else
        loss_value_old = loss_value_new;
    end
    
end

score = loss_value_new;



function f_val = function_eval(W_obs, W_dir, lambda)
% Calculate the value of objective function (3)
W_obs_hat = W_dir+ W_dir*W_dir;
W_obs_hat = W_obs_hat - diag(diag(W_obs_hat));
f_val =  sum(sum((W_obs-W_obs_hat).^2)) + lambda*sum(sum(W_dir.^2));