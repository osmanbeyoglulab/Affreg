%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figure 4c, S16 RBP contacting residue importance maps
%

function run_rbp_projection

model.lambda = 1e-3;
model.spectrumA = 0.9;
model.spectrumB = 1;
model.rsL2 = 0;

load('data/RbpData_1.0_norep.mat');
dataset = RbpData;
ss = padstr(dataset.Proteins.seq);
names = RbpData.Proteins.name;
seqs = dataset.Proteins.seq;


% train a classifier with strings split
% break strings if they contain domains
scat = [];
W = []; N=[];
k = 1;IX = [];
for i=1:length(seqs)
    s = seqs{i};
    ix = [0 strfind(s,'*') length(s)+1];
    for j=1:length(ix)-1
        scat{k} = s(ix(j)+1:ix(j+1)-1);
        N{k} = sprintf('%s.domain_%d', names{i}, j);
        k = k + 1;
    end
    IX = [IX; repmat(i, length(ix)-1, 1)];
end
sscat = padstr(scat);

% structure dataset to contain multiple domains
Y_train = quantilenorm(RbpData.intensities');
Y = bsxfun(@times, Y_train, 1./sqrt(sum(Y_train.^2)));
Y1 = Y(:, IX);

D = RbpData.kmer_counts;
D1 = bsxfun(@minus, RbpData.kmer_counts, mean(D));
D_norm = bsxfun(@times, D1, sqrt(1./sum(D1.^2)));

for i=1:size(sscat, 1); 
    P1(i, :) = (~cellfun(@isempty, regexp(sscat(i, :), dataset.Proteins.kmers)))'; 
end;
P_norm = bsxfun(@times, P1, 1./sqrt(sum((P1+eps).^2)));

fullmodel = ar_train(D_norm,P_norm,Y1,0.1, model.rsL2, model.spectrumA, 1);

w = ar_model2w(fullmodel);
kmer_projections1 = Y1'*D_norm * w;
for six = 1:size(sscat,1)
    [sw1(six, :), w_ix{six}] = proj2seqw( sscat(six, :), kmer_projections1(six, :), dataset.Proteins.kmers);
end


% plot the rbp selected heatmap
ix_snrpa1 =find(~cellfun(@isempty, strfind(N, 'SNRPA'))); 
ix_hur1 =find(~cellfun(@isempty, strfind(N, 'HuR')));
ix_pab1 =find(~cellfun(@isempty, strfind(N, 'PABPC1')));
ix_rbfox1 =find(~cellfun(@isempty, strfind(N, 'RBFOX1')));
ix_fox1 =find(~cellfun(@isempty, strfind(N, 'FOX-1')));
ix_paball1 =find(~cellfun(@isempty, strfind(N, 'PABP')));
ix_foxall1 =find(~cellfun(@isempty, strfind(N, 'FOX')));
ix_igfall1 =find(~cellfun(@isempty, strfind(N, 'IGF')));
select = [ix_foxall1, ix_snrpa1, ix_paball1, ix_hur1, ix_igfall1];

for six = 1:size(sscat,1)
    [swselect(six, :)] = proj2seqw( sscat(six, :), kmer_projections1(six, :), dataset.Proteins.kmers, [], 1);
end
h=figure; plotGrid(swselect(select,:), sscat(select, :), N(select), 10);
colormap jet
c=colormap();
colorbase = c(round(4:1.5:end),:); 
c1=[colorbase; repmat(c(end, :), 64-size(colorbase, 1), 1)];
colormap(c1)
colorbar
set(gca, 'xticklabel',{})
set(h, 'PaperUnits', 'inches', 'PaperSize',[20 50], 'PaperPosition', [0 0 20 50])
%set(h, 'Units', 'pixels', 'Position', [0, 0, 2500, 1000], 'PaperPositionMode','auto','visible','off')
%print(h, sprintf('results/map_rbpselect_sw1_splitdomain.png'),'-dpng', '-r300')


% plot rbp heatmap
h=figure; plotGrid(swselect, sscat, N, 7);
colormap jet
c=colormap();
colorbase = c(round(6:1.7:end),:); 
c1=[colorbase; repmat(c(end, :), 64-size(colorbase, 1), 1)];
colormap(c1)
colorbar
set(gca, 'xticklabel',{})
%set(h, 'Units', 'pixels', 'Position', [0, 0, 3500, 2000], 'PaperPositionMode','auto','visible','off')
%set(h, 'PaperUnits', 'inches', 'PaperSize',[20 50], 'PaperPosition', [0 0 20 50])
%print(h, sprintf('results/map_rbp_sw1_splitdomain.png'),'-dpng', '-r300')


function pad = padstr(s)
maxlen = max(cellfun(@length, s));
pad = char(zeros(length(s), maxlen) + '-');
for i=1:length(s)
    pad(i, 1:length(s{i})) = s{i};
end


function plotProjection(seq, w, Rnd)
l = length(seq);
if ischar(seq)
    seq = cellstr(seq(:));
end
hold on
plot(Rnd','r');
plot(w,'k', 'linewidth',2);
set(gca, 'xtick',1:l); 
set(gca, 'xticklabel',seq);
