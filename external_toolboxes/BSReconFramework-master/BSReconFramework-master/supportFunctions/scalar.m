function N = scalar(A)
    N = sum(sum(A .* conj(A),2),1);
end