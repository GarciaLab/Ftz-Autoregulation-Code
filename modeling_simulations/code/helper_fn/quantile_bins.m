function [bins,dist_quant_edges] = quantile_bins(X,N)
N = N-1;
dist_quant_edges = [min(X),quantile(X,N),max(X)];
bins = discretize(X,dist_quant_edges);
end