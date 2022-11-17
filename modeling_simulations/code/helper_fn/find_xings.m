function cross_ix = find_xings(x,y)

x = x(:);
y = y(:);

if numel(x) ~= numel(y)
    error('x and y must be the same size');
end

cross_ix = find(diff(sign(x - y)) ~= 0);

% remove neighboring indices
dix = diff(cross_ix);
cross_ix([dix == 1;false]) = [];

end