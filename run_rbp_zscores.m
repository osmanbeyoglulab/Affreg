%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figures 4e, S19, S20 predictions for RBP motifs (zscore model)
%

function run_rbp_zscores
source = 'rbpz';

is_centered = 1;
isexpy = 0;
fprintf('is centered = %d, isexpy = %d\n', is_centered, isexpy)
fprintf('parameters set to win top 100 precision, and therefore produce better motifs') 

[names, seqs] = fastaread('data/rbpz_protseq.fa');
D = csvread('data/rbpz_probe6_counts.csv',1,1);
P = csvread('data/rbpz_prot_counts.csv',1,1);
Y = dlmread('data/rbpz_Y.txt','\t', 1,1);
Y = quantilenorm(zscore(Y));

D = double(D>0);
P = double(P>0);
P = P(:,sum(P)>2); % keep kmers that occur more than twice
Dn = bsxfun(@times, D, 1./sqrt(sum(D.^2))); %unit norm D
D = Dn;
Y_gr = Y;
expfactor = 1.1;

if isexpy
    Y = (expfactor).^Y-1; % exponentiate and center (exp(0)-1=0)
end

if is_centered
    mu = mean(Y,2);
else
    mu = zeros(size(Y, 1), 1);
end
Y = Y - repmat(mu, 1, size(Y,2));

n = size(Y,2); %number of proteins

n_folds = 10;
fold_ix = round(linspace(0,n, n_folds+1));

% set lambdas, rsL2 and spectrum as vectors for parameter search
lambdas = [1]; 
rsL2 = 0;
spectrumA = [0.9];
spectrumB = [0.9];
[L1 L2 SA SB] = ndgrid(lambdas, rsL2, spectrumA, spectrumB);
params = [L1(:) L2(:) SA(:) SB(:)];

%parfor param_ix = 1:size(params,1)
base = 0;

models = cell(n_folds, size(params, 1));
preds = cell(n_folds, size(params, 1));
for param_ix = 1:size(params,1)
    
    lambda = params(param_ix, 1);
    rsL2 = params(param_ix, 2);
    spectrumA = params(param_ix, 3);
    spectrumB = params(param_ix, 4);
    fprintf('%d: starting %d/%d\n',param_ix, param_ix, size(params,1));
    
    models_par = cell(n_folds, 1);
    for fold = 1:n_folds
        %fprintf('fold %d/%d',fold, n_folds);
        fprintf('%d: specA=%.3f specB=%.3f lambda=%.3f rsL2=%.3f fold=%d\n', param_ix, spectrumA, spectrumB, lambda, rsL2, fold);
        
        % create subsets
        test_ix = fold_ix(fold)+1:fold_ix(fold+1);
        train_ix = setdiff(1:n, test_ix);
        
        % create subsets
        Y_test  = Y(:, test_ix);
        Y_train = Y(:, train_ix);
        P_test  = P(test_ix, :);
        P_train = P(train_ix,:);
        seqs_test = seqs(test_ix);
        seqs_train = seqs(train_ix);
        
        models{fold, param_ix} = ar_train(D, P_train, Y_train, lambda, rsL2, spectrumA, spectrumB);
        preds{fold, param_ix}  = ar_predict(D, P_test, Y_train, models{fold, param_ix});

        preds_nn{fold, param_ix} = predictNN(P_test, P_train, Y_train);
        preds_or{fold, param_ix} = predictBN(Y_test, Y_train);
        preds_bl{fold, param_ix} = predictBL(seqs_test, seqs_train, Y_train);
    end
end

ar_rec = zeros(size(Y, 1), size(Y, 2), size(params, 1));
nn_rec = zeros(size(Y, 1), size(Y, 2), size(params, 1));
bl_rec = zeros(size(Y, 1), size(Y, 2), size(params, 1));
or_rec = zeros(size(Y, 1), size(Y, 2), size(params, 1));
for param_ix = 1:size(params,1)
    fprintf('%d: reconstructing %d/%d\n',param_ix, param_ix, size(params,1));
    for fold = 1:n_folds
        % create subsets
        test_ix = fold_ix(fold)+1:fold_ix(fold+1);
        ar_rec(:, test_ix, param_ix) = preds{fold, param_ix}.rec + repmat(mu, 1, size(test_ix, 2));
        nn_rec(:, test_ix, param_ix) = preds_nn{fold, param_ix} + repmat(mu, 1, size(test_ix, 2));
        bl_rec(:, test_ix, param_ix) = preds_bl{fold, param_ix} + repmat(mu, 1, size(test_ix, 2));
        or_rec(:, test_ix, param_ix) = preds_or{fold, param_ix} + repmat(mu, 1, size(test_ix, 2));
    end
end

% max_ix picks the best models
% can be set for example as [m max_ix] = max(diag(corr(Y_gr, ar_rec)))
% if multiple parameter settings are tested
max_ix = 1;

folder = performancePlots(ar_rec(:,:, max_ix), nn_rec(:,:, max_ix), bl_rec(:,:, max_ix), or_rec(:,:, max_ix), Y_gr, source, 'run_rbp_zscores.m', 0, 'spearman');
disp('run_rbpz_zscores.m done.')

