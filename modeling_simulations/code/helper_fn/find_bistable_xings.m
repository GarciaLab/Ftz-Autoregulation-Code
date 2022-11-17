function [x,f] = find_bistable_xings(fun,guesses)

guesses(isnan(guesses) | isinf(guesses)) = [];

x = nan(size(guesses));
f = nan(size(guesses));

if numel(guesses) <= 3
    for ii = 1:length(guesses)
        if ~isnan(fun(guesses(ii))) && isreal(fun(guesses(ii))) && ...
                ~isinf(fun(guesses(ii)))
            [x(ii),f(ii),~,~] = fzero(fun,guesses(ii));
        end
    end
end

end