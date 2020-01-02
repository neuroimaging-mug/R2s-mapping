function [pattern,R_True] = gen_random_pattern(pdf,R,wind)

% [PATTERN, DENSITY, R_TRUE, PDF] = GEN_DATA_RANDOM_TEMPLATE(TEMPLATE,R)
%
% INPUT:
% pdf:      sampling pdf
% R:        desired reduction factor
% wind:     window size for filling of holes
%
% OUTPUT
% pattern:   sampling pattern
% R_True:    true reduction factor

%% generate optimized random subsampling pattern
[N,M] = size(pdf);

% normalize pdf (int(pdf)=1)
pdf = pdf/sum(pdf(:));
P   = pdf(:);

% cumulative distribution
P   = cumsum(P) ./ sum(P);
S   = floor(length(P)/R); % safety margin since indices can be redrawn

% try several random patterns, take most incoherent one
pattern = zeros(size(pdf));
minval  = 1e99;
maxit   = 6;

for n=1:maxit
    mask  = [];
    draws = S;
    while draws>0
        d_mask = zeros(draws,1);
        for i = 1:draws
            X = rand < P;
            d_mask(i) = sum(X);
        end
        d_mask = (length(P) + 1) - d_mask;
        mask   = unique([mask;d_mask]);
        draws  = S-length(mask);
    end
    
    try_pattern = zeros(size(pdf));
    try_pattern(mask) = 1;
    
    coherence = ifft2(try_pattern./pdf);
    max_coher = max(abs(coherence(2:end)));
    if  max_coher < minval
        minval  = max_coher;
        pattern = try_pattern;
    end
end

% fill holes
for i=1:N
    for j = 1:M
        imin = max(1,i-wind); imax = min(N,i+wind);
        jmin = max(1,j-wind); jmax = min(M,j+wind);
        if max(pattern(imin:imax,jmin:jmax)) < 1 % window empty
            pattern(i,j) = 1;
        end
    end
end

% effective acceleration factor
R_True = length(P)/nnz(pattern);

% end main function


