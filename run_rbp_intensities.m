%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figure 4a,b and S14 RBP intensity performance plots
%

function run_rbp_intensities

version = 'rbp_intensities';

stream = RandStream.getGlobalStream;
reset(stream,100);

load('data/RbpData_1.0_norep.mat');
n = size(RbpData.intensities,1); %number of proteins
r = randperm(n);

% permute the data
RbpData.motifnames = RbpData.motifnames(r);
RbpData.intensities = RbpData.intensities(r, :);
RbpData.Proteins.kmer_counts = RbpData.Proteins.kmer_counts(r,:); 
RbpData.Proteins.uniqueid = RbpData.Proteins.uniqueid(r);
RbpData.Proteins.motifid = RbpData.Proteins.motifid(r);
RbpData.Proteins.geneid = RbpData.Proteins.geneid(r); 
RbpData.Proteins.name = RbpData.Proteins.name(r); 
RbpData.Proteins.replicateid = RbpData.Proteins.replicateid(r);
RbpData.Proteins.seq = RbpData.Proteins.seq(r);

% set the basic data elements
Y = RbpData.intensities';
P_norm = RbpData.Proteins.kmer_counts;
D = RbpData.kmer_counts;
D = bsxfun(@minus, D, mean(D));  % centering D improves oracle. Not necessary for AR to win.
D_norm = bsxfun(@times, D, sqrt(1./sum(D.^2)));
seq = RbpData.Proteins.seq;
ids = RbpData.Proteins.uniqueid;

n = size(Y,2); %number of proteins
n_folds = 10;
fold_ix = round(linspace(0,n, n_folds+1));

lambda = 1.0000e-03;
rsL2 = 0;
spectrumA = 1;
spectrumB = 0.85;

Pred_rec = zeros(size(Y));
Pred_nn = zeros(size(Y));
Pred_bl = zeros(size(Y));
Pred_bn = zeros(size(Y));
Pred_test = zeros(size(Y));

for fold = 1:n_folds
    test_ix = fold_ix(fold)+1:fold_ix(fold+1);
    train_ix = setdiff(1:n, test_ix);
    
    % create subsets
    Y_test = Y(:,test_ix);
    Y_train = Y(:, train_ix);
    
    %q-normalize with Y_train only, inject values to Y_test
    Y_train = quantilenorm(Y_train);
    Y_train = bsxfun(@times, Y_train, 1./sqrt(sum(Y_train.^2)));
    norm_values = sort(Y_train(:,1));
    if 1
        for i = 1:size(Y_test,2)
            [~, sort_ix] = sort(Y_test(:,i));
            Y_test(sort_ix,i) = norm_values;
        end
    end
    
    P_test = P_norm(test_ix,:);
    P_train = P_norm(train_ix,:);
    
    mu = mean(Y_train, 2);
    mu = zeros(length(mu), 1); %no need to center for original results, although centering can improve performance
    Y_center = bsxfun(@minus, Y_train, mu);
    
    model = ar_train(D_norm,P_train,Y_center,lambda, rsL2, spectrumA, spectrumB, 1);
    pred = ar_predict(D_norm, P_test, Y_center, model);

    % predict
    %   nn - nearest neighbor
    %   bn - best neighbor (oracle)
    %   bl - blosum neighbor
    P_test_seq = RbpData.Proteins.seq(test_ix);
    P_train_seq = RbpData.Proteins.seq(train_ix);
    [pred_nn, ~] = predictNN(P_test, P_train, Y_train);
    [pred_bn, ~] = predictBN(Y_test, Y_train);
    [pred_bl, ~] = predictBL(P_test_seq, P_train_seq, Y_train);
    
    %save predicted probe intensities
    Pred_rec(:,test_ix) = pred.rec + repmat(mu, 1, size(pred.rec, 2));
    Pred_nn(:,test_ix) = pred_nn;
    Pred_bl(:,test_ix) = pred_bl;
    Pred_test(:,test_ix) = Y_test;
    Pred_bn(:,test_ix) = pred_bn;
    
    fprintf('finished fold %d/%d\n', fold, n_folds);
end
%performancePlots(Pred_rec, Pred_nn, Pred_bl, Pred_bn, Pred_test, version, 'run_rbp_intensities.m', 0, 'pearson');
1;


