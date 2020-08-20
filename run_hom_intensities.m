%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figure 1f,g performance plots
%

% function run_hom_intensities
close all
clear variables

%fix the random seed for figure 
stream = RandStream.getGlobalStream;
reset(stream,100);

% % load and setup data
% load("C:/Users/XIM33/Documents/Dr. Osmanbeyoglu/Project/affreg/PythonVersion/github/data/10X_pbmc5k_seurat_omni1.mat");
% D = double(D); %gene X TF
% P = Ppbmc5kcd14; %protein X sample
% Y = Ypbmc5kcd14;% gene X sample


load 'data/PbmDataHom6_norm.mat';
P = PbmData.Proteins.kmer_counts;
% % sample p columns
% nsamples2 = size(P,2);
% samples2 = randi([1 nsamples2],1,round(nsamples2/2));
% P = P(:,samples2);

%% sample P rows
% nsamples1 = size(P,1);
% samples1 = randi([1 nsamples1],1,round(nsamples1/2));
% P = P(samples1,:)
%sample end

D = PbmData.kmer_counts;
Y = PbmData.intensities;
% sample Y
% Y = Y(:, samples1);

% % randomize D
D1 = size(D,1) ; % example matrix
D2 = size(D,2) ;
for k = 1:D2
   D(:,k) = D(randperm(D1),k) ;
end
% for k = 1:D1
%    D(k,:) = D(k,randperm(D2)) ;
% end
n = size(Y,2); %number of proteins
r = randperm(n);
P = P(r,:);
Y = Y(:,r);
% P_seq = PbmData.Proteins.seqs(r);

% % randomize P
% P1 = size(P,1) ; % example matrix
% P2 = size(P,2) ;
% % for k = 1:P2
% %    P(:,k) = P(randperm(P1),k) ;
% % end
% for k = 1:P1
%    P(k,:) = P(k,randperm(P2)) ;
% end



%names = PbmData.names(r);
P_norm = P;
D_norm = D;

% setup for 10 fold cross validation
n_folds = 10;
fold_ix = round(linspace(0,n, n_folds+1));

% containers for learning predictions
Pred_rec = zeros(size(Y));
Pred_nn = zeros(size(Y));
Pred_bn = zeros(size(Y));
Pred_bl = zeros(size(Y));
Pred_test = zeros(size(Y));

% optimal learning setting
lambda = 0.01;
rsL2 = 0;
spectrumA = 1;
spectrumB = 0.8;

for fold = 1:n_folds
    % create subsets
    test_ix = fold_ix(fold)+1:fold_ix(fold+1);
    train_ix = setdiff(1:n, test_ix);
    
    % create subsets
    Y_test = Y(:,test_ix);
    %names_test = names(:,test_ix);
    Y_train = Y(:, train_ix);
    %names_train = names(train_ix);
    
    %normalize Y_train, and inject values onto Y_test, keeping them
    %separate
    Y_train = quantilenorm(Y_train);
    Y_train = bsxfun(@times, Y_train, 1./sqrt(sum(Y_train.^2)));
    norm_values = sort(Y_train(:,1));
    if 1
        for i = 1:size(Y_test,2)
            [~, sort_ix] = sort(Y_test(:,i));
            Y_test(sort_ix,i) = norm_values;
        end
    end
%     Y_test = quantilenorm(Y_test);
%     Y_test = bsxfun(@times, Y_test, 1./sqrt(sum(Y_test.^2)));    

    % exponentially scale the distribution to increase the importance of
    % top probes
    Y_test = exp(100*(Y_test))-1;
    Y_test = bsxfun(@times, Y_test, 1./sqrt(sum(Y_test.^2)));
    Y_train = exp(100*(Y_train))-1;
    Y_train = bsxfun(@times, Y_train, 1./sqrt(sum(Y_train.^2)));
    
    P_test = P_norm(test_ix,:);
    P_train = P_norm(train_ix,:);

     
    % affinity regression train and predict
    model = ar_train(D_norm,P_train,Y_train,lambda, rsL2, spectrumA, spectrumB);
    new_pred5 = ar_predict(D_norm, P_test, Y_train, model);

    % predict using the competitors
    % pred_nn: nearest neighbor
    % pred_bn: best nearest neighbor (oracle)
    % pred_bl: blosum nearest neighbor
    
%     P_test_seq = P_seq(test_ix);
%     P_train_seq = P_seq(train_ix);
    [pred_nn, ~] = predictNN(P_test, P_train, Y_train);
    [pred_bn, ~] = predictBN(Y_test, Y_train);
%     [pred_bl, ~] = predictBL(P_test_seq, P_train_seq, Y_train);
       
    %save predicted probe intensities
    Pred_rec(:,test_ix) = new_pred5.rec;
    Pred_nn(:,test_ix) = pred_nn;
    Pred_bn(:,test_ix) = pred_bn;
%     Pred_bl(:,test_ix) = pred_bl;
    Pred_bl(:,test_ix) = pred_bn;
    Pred_test(:,test_ix) = Y_test;
    
    fprintf('fold %d/%d\n',fold, n_folds);
end
% print(corr(Pred_rec(:),Pred_test(:)));
performancePlots(Pred_rec, Pred_nn, Pred_bl, Pred_bn, Pred_test, 'hom_intensities', 'run_hom_intensities.m', 1, 'pearson');

