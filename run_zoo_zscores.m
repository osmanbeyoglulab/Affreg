%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figure 3d train on alleyene and predict on zoo performance plots (zscore model)
%

function run_zoo_zscores

% learn a zscore model from alleyene data, and predict on zoo zscores
source = 'zoo218';

% load zoo data
zoo_counts = importdata('data/zoo_counts.csv');
zoo_zscores = importdata('data/zoo_y_nodup.csv');
P_zoo = zoo_counts.data;
y_zoo = zoo_zscores.data;
seqs = importdata('data/zoo_seqs_218.fa');
%names_zoo = seqs(1:2:end);
seqs_zoo  = seqs(2:2:end);

% load alleyene data for training
data = load('data/alleyene_data.mat');
prot_counts = importdata('data/alleyene_prot_counts.csv');
seqs_al = importdata('data/alleyene_prots_seqs.fa');
%names_al = seqs_al(1:2:end);
seqs_al  = seqs_al(2:2:end);

P = double(prot_counts.data>0);
ix = sum(P)>3;
P = P(:, ix);
P_zoo = double(P_zoo(:, ix) > 0);
D = double(data.D>0);
ix = sum(D)>3;
D = D(:, ix);

Y_al = quantilenorm(zscore(data.Y));
Y_algr = data.Y;
Y_zoo = y_zoo;

isexpy = 0;
is_centered = 1;
expfactor = 1.5;

if isexpy
    Y_al = expfactor.^(Y_al);
end

if is_centered
    mu = mean(Y_al,2);
else
    mu = zeros(size(Y_al, 1), 1);
end
Y_al = Y_al - repmat(mu, 1, size(Y_al,2));


n = size(Y_al,2); %number of proteins

lambda = 1;
rsL2 = 0;
spectrumA = 1;
spectrumB = 1;

fprintf('starting %f\n', lambda);
model = ar_train(D, P, Y_al, lambda, rsL2, spectrumA, spectrumB);
zoo_preds = ar_predict(D, P_zoo, Y_al, model);
%preds = ar_predict(D, P, Y_al, model);  % evaluate learned model

% get predictions for
% nn - nearest neighbor
% or - oracle neighbor
% bl - blosum neighbor
idx_nn = knnsearch(P, P_zoo);
idx_or = knnsearch(Y_algr', Y_zoo');
similarity = nwalign(seqs_zoo, seqs_al);
[~, idx_bl] = max(similarity,[],2);

nn_rec = Y_algr(:,idx_nn); % + repmat(mu, 1, length(idx_nn));
or_rec = Y_algr(:,idx_or); % + repmat(mu, 1, length(idx_nn));
bl_rec = Y_algr(:,idx_bl);% + repmat(mu, 1, length(idx_nn));

% get the affinity regression prediction
ar_rec = zoo_preds.rec + repmat(mu, 1, length(idx_nn));

performancePlots(ar_rec, nn_rec, bl_rec, or_rec, Y_zoo, source, 'run_zoo_zscores.m', 1,'spearman');

disp('run_zoo_zscores done.')

