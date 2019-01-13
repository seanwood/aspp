function [H] = inferCodebookCoefficients(V, W, maxUpdates, alpha, epsilon, seedValue)

rand('seed', seedValue);
    
dictionarySize = size(W, 2);
numExamples = size(V, 2);
H = rand( dictionarySize, numExamples ) + epsilon;

for updateIndex = 1:maxUpdates
    H = H .* ( ( W.' * ( V ./ (W * H) ) ) ./ ( sum(W, 1).' + alpha + epsilon ) );
end

end

