function [H_star, H_c_star, H_s_star, final_H_c, final_H_s] = PSMVC(Data_set, K, lambda, eta, repeat_times, iter, rho)

% Estimating model parameters Hc and Hs by updating them
% according to Equations (4), and (5) alternately. It repeats the entire
% updating procedure multiple times and chooses the result that gives
% the lowest value of objective function in Equation (3) as the final estimator.


% Input:
%   Data_set: a cell, in which each input data is named by
%   Data_set{l}, and the score matrix for this data is stored in Data_set{l}.network; 
%   K: maximum number of possible protein complexes. The default value is
%   2500;
%   lambda: parameter which controls the effect of low rank constrain. The
%   default value is 2^2;
%   eta: the common factor ratio. The default value is 0.5;
%   repeat_times: the number of times that we repeat the entire calculation
%   to avoid local minimum. The default value is 20;
%   iter: the number of iterations limited in PSMVC. The default value is
%   200;
%   rho: the tolerance threshold of the stop criterion. The default value is 1e-6.




% Outputs:
%   H_star: the overall protein-complex assignment matrices.
%   H_c_star: the shared complexes between different data sources.
%   H_s_star: the specific complexes within each data source.
%   final_H_c: the estimator of Hc.
%   final_H_s: the estimator of Hs

if nargin < 7
    rho = 1e-6;
end

if nargin < 6
    iter = 200;
end

if nargin < 5
    repeat_times = 20;
end

if nargin < 4
    eta = 0.5;
end


if nargin <3
    lambda = 4;
end

if nargin < 2
    K = 2500;
end

if nargin < 1
    error('You need input dataset');
end

fprintf('Estimating parameters...')
fprintf('\n')

 
% Calculate the number of data sources.
L = length(Data_set);

% Calculate the dimension of common latent factor and specific latent factor.
Kc = eta*K;
Ks = (1/L)*(K - Kc);

% Calculate the number of proteins.
N = size(Data_set{1}.network,1);

lowest_score = inf;
final_score = inf;
% Repeat the entire updating procedure multiple times.
for i = 1: repeat_times
    fprintf(['This is the ',num2str(i), '-th repeat...'])
    fprintf('\n')

% Initialize matrix Hc    
    H_c_old = rand(N,Kc);
    
% Initialize matrix Hs and calculate theta
    for l = 1 : L
        Data_set{l}.network = Data_set{l}.network - diag(diag(Data_set{l}.network));
        H_s_old{l} = rand(size(Data_set{l}.network,1),Ks);
        
% Calculate theta        
        theta{l} =zeros(size(Data_set{l}.network,1),1);
        ID = sum(Data_set{l}.network) > 0;
        theta{l}(ID) = 1;
        theta{l} = theta{l}*theta{l}';
        
    end
    score = zeros(1,iter);
% Update Hc and Hs according to the updating rules.
    for j  = 1: iter

        Ws = cell(L,1);
        Ws(:) = {0};

        Wc = ((theta{1}.*Data_set{1}.network)./([H_c_old H_s_old{1}]*[H_c_old H_s_old{1}]' + eps))*H_c_old + ((theta{2}.*Data_set{2}.network)./([H_c_old H_s_old{2}]*[H_c_old H_s_old{2}]' + eps))*H_c_old;

        for l = 1 : L
            Ws{l} = ((theta{l}.*Data_set{l}.network)./([H_c_old H_s_old{l}]*[H_c_old H_s_old{l}]' + eps))*H_s_old{l};
            H_s{l} = 0.5*H_s_old{l} + 0.5*H_s_old{l}.*((Ws{l})./(theta{l}*H_s_old{l} + lambda*H_s_old{l} + eps));
        end
        
        H_c = 0.5*H_c_old + 0.5*H_c_old.*((Wc)./(theta{1}*H_c_old + theta{2}*H_c_old + lambda*H_c_old + eps));

       score(j) = - sum(sum(theta{1}.*Data_set{1}.network.*log([H_c H_s{1}]*[H_c H_s{1}]' + eps))) + sum(sum(theta{1}.*([H_c H_s{1}]*[H_c H_s{1}]')))  - sum(sum(theta{2}.*Data_set{2}.network.*log([H_c H_s{2}]*[H_c H_s{2}]' + eps))) + sum(sum(theta{2}.*([H_c H_s{2}]*[H_c H_s{2}]'))) + lambda*sum(sum(H_c.*H_c)) + lambda*sum(sum( H_s{1}.*H_s{1})) + lambda*sum(sum( H_s{2}.*H_s{2}));
        
        if abs(score(j) - lowest_score)/abs(lowest_score) < rho
            break;
        else        
            H_c_old = H_c;
            H_s_old = H_s;
            lowest_score = score(j);
        end
    end
    
        if score(j) < final_score
        final_H_c = H_c;
        final_H_s = H_s;
        final_score = score(j);
        end
end

% Detect protein complexes from the estimators of model parameters.

fprintf('Detecting protein complex ...')
fprintf('\n')
H_star = final_H_c;

% Obtain H_star
for l = 1 : L
    H_star = [H_star, final_H_s{l}];
end

save score2 score
    H_star(H_star < 10^(-3)) = 0;
    [H_sort,~] = sort(H_star,2,'descend');
    dif_H = H_sort(:,1:(size(H_sort,2) - 1)) - H_sort(:,2:size(H_sort,2));
    [~,dif_ID] = max(dif_H,[],2);
    act_val = zeros(N,1);
    for i = 1 : N
       act_val(i) = H_sort(i,dif_ID(i) + 1);
    end
    act_matrix = repmat(act_val,1,size(H_star,2));
    H_act = zeros(N,size(H_star,2));
    H_act(H_star > act_matrix) = 1;
    H_star = H_act;
    H_c_star = H_star(:,1:Kc);
    H_s_star = H_star(:,(Kc+1):size(H_star,2));
    
    % Filter out the detected complexes which include less than three proteins.
    small_size_indices = sum(H_c_star) <= 2;
    H_c_star = H_c_star(:,~small_size_indices);
    small_size_indices = sum(H_s_star) <= 2;
    H_s_star = H_s_star(:,~small_size_indices);
    H_star = [H_c_star H_s_star];
   