%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%
%       generate figure 4c, S16 RBP contacting residue importance maps
%

function run_hom_projection
clrred = [228,26,28]/255;
clrblue = [55,126,184]/255;
clrgreen = [77,175,74]/255;
clrpurple = [152,78,163]/255;
clrblack = [0 0 0 ];

model.lambda = 0.001;
model.rsL2 = 0;
model.spectrumA = 1;
model.spectrumB = 0.8;


load 'data/PbmDataHom6_norm.mat';
P = PbmData.Proteins.kmer_counts;
D = PbmData.kmer_counts;
Y = PbmData.intensities;
Y_norm = quantilenorm(Y);
Y_train = bsxfun(@times, Y_norm, 1./sqrt(sum(Y_norm.^2)));
Y_skew100 = exp(100*(Y_train))-1; 
Y_skew100 = bsxfun(@times, Y_skew100, 1./sqrt(sum(Y_skew100.^2)));
    
Y = Y_train;
Y = Y_skew100;

ix_pkn      = find(~cellfun(@isempty, strfind(PbmData.names,'Pknox1')));
ix_hoxa9    = find(~cellfun(@isempty, strfind(PbmData.names,'Hoxa9')));
ix_pyp      = find(~cellfun(@isempty, strfind(PbmData.Proteins.seqs,'PYP')));
ix_hox      = find(~cellfun(@isempty, strfind(PbmData.names,'Hox')));
ix_hox      = [ix_hox(:); ix_pyp(:)];

fullmodel = 	(D,P,Y,model.lambda, model.rsL2, model.spectrumA, model.spectrumB, 1);
w = ar_model2w(fullmodel);
kmer_projections = Y' * D * w;

%map projection onto sequences
seqs =  cellstr(multialign(PbmData.Proteins.seqs));
chrs = char(seqs);
dashs = chrs == '-';
sds_ix = sum(dashs)<10;
full_chrs = chrs(:, sds_ix);

W  = proj2seqw(full_chrs, kmer_projections, PbmData.Proteins.kmers);
Whox = proj2seqw(full_chrs(ix_hox,:), kmer_projections(ix_hox,:), PbmData.Proteins.kmers);
Wpyp = proj2seqw(full_chrs(ix_pyp,:), kmer_projections(ix_pyp,:), PbmData.Proteins.kmers);

% plot projection maps
figure, plotGrid(W, full_chrs, PbmData.names); colormap jet
figure, plotGrid(Whox, full_chrs(ix_hox, :), PbmData.names(ix_hox)); colormap jet
figure, plotGrid(Wpyp, full_chrs(ix_pyp, :), PbmData.names(ix_pyp)); colormap jet

% plot projection plots for pknox1 and hox
h = figure; hold on
plot ( W(ix_pkn,:),'color', clrblue); title ('Positional importance')
plot ( W(ix_hoxa9,:),'color',clrred);

for i=1:size(W, 2); labs{i} = full_chrs(ix_pkn, i); end;
set(gca, 'xtick', 1:size(W, 2)); %set(gca, 'xticklabel', labs );
h_legend = legend('Pknox1','HoxA9','location','northwest');
set(h_legend,'FontSize',14);
for i=1:size(W, 2); labs_hoxa9{i} = full_chrs(ix_hoxa9, i); end;
%set(gca, 'xtick', 1:size(W, 2)); set(gca, 'xticklabel', labs_hoxa9 );
set(gca, 'xticklabel', [])
xtick = 0;
ytick =0;
for i =1:size(W, 2); text(xtick+i-1.2, ytick-0.002, labs{i},'color',clrblue, 'fontname','courier'); end
for i =1:size(W, 2); text(xtick+i-1.2, ytick-0.005, labs_hoxa9{i},'color',clrred, 'fontname','courier'); end

%plot the helixes
yoff = -0.007;
line([1 size(W, 2)+1],[yoff yoff],'color',clrblack)
rectangle('Position',[8.5,-0.001+yoff,11.5,0.002],'curvature', 0.1,'FaceColor', clrpurple)
rectangle('Position',[26.5,-0.001+yoff,10,0.002], 'curvature', 0.1,'FaceColor', clrpurple)
rectangle('Position',[41.5,-0.001+yoff,15,0.002], 'curvature', 0.1,'FaceColor', clrpurple)
text(14, yoff+0.003,'\alpha1')
text(31, yoff+0.003,'\alpha2')
text(49, yoff+0.003,'\alpha3')

set(gca, 'xtick', [])
%set(gca, 'ytick', [])
ylabel('Importance')

[wmat handle] = seqlogo(cellstr(full_chrs), 'alphabet', 'NT');
[wmat handle] = seqlogo(cellstr(full_chrs(ix_hox,:)), 'alphabet', 'NT');
[wmat handle] = seqlogo(cellstr(full_chrs(ix_pyp,:)), 'alphabet', 'NT');

1;


function seqW = proj2seqw (seqs, proj, kmers)
seqW = zeros(size(seqs));
n_kmer = length(kmers{1});
for i=1:size(seqs,1)
    seq = seqs(i,:);
    seq_ix = bsxfun(@plus, repmat(0:(n_kmer-1), length(seq)-n_kmer+1, 1), (1:length(seq)-n_kmer+1)');
    %seq_ix = seq_ix + 1;
    seq_kmers = cellstr(seq(seq_ix));
    w_ix = cellfun(@(x) find(strcmpi(kmers, x)), seq_kmers,'uniformoutput',0);
    empty_ix = cellfun(@isempty, w_ix);
    seqW(i, ~empty_ix) = cell2mat(cellfun(@(x) proj(i,x), w_ix','uniformoutput',0));

    w = seqW(i,:);
    w1 = [w(1) w(1:end-1)];
    w2 = [w(1:2) w(1:end-2)];
    w3 = [w(1:3) w(1:end-3)];

    seqW(i,:) = (w+w1+w2+w3)/4;
end



