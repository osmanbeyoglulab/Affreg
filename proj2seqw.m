%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rafi Pelossof, MSKCC, 2015
%

function [seqW, w_ix] = proj2seqw (seqs, proj, kmers, w_ix, isavg)
% seqs: char array, each row is a sequence
% proj: projections, #rows = #rows(seqs), #cols=#rows(kmers)
% kmers: cellstr vector, #row=all kmers, #col=1
if nargin < 4
    w_ix = [];
    isavg = false;
end

if nargin < 5
    isavg = false;
end


seqW = zeros(size(seqs));
n_kmer = length(kmers{1});
for i=1:size(seqs,1)
    seq = seqs(i,:);
    seq_ix = bsxfun(@plus, repmat(0:(n_kmer-1), length(seq)-n_kmer+1, 1), (1:length(seq)-n_kmer+1)');
    %seq_ix = seq_ix + 1;
    seq_kmers = cellstr(seq(seq_ix));
    if (isempty(w_ix))
        w_ix = cellfun(@(x) find(strcmpi(kmers, x)), seq_kmers,'uniformoutput',0);
    end
    empty_ix = cellfun(@isempty, w_ix);
    seqW(i, ~empty_ix) = cell2mat(cellfun(@(x) proj(i,x), w_ix','uniformoutput',0));
    
    w = seqW(i,:);
    w1 = [w(1) w(1:end-1)];
    w2 = [w(1:2) w(1:end-2)];
    w3 = [w(1:3) w(1:end-3)];
    
    if isavg
        seqW(i,:) =  (w+w1+w2+w3)/4;
    else
        seqW(i,:) = w;
    end
end
