function j=j_sample(n, jamp)

% Return an asymmetric interaction matrix with random Gaussian entries
% and diagonal entries set to zero

j=jamp/sqrt(n)*randn(n);
j=j-diag(diag(j));
