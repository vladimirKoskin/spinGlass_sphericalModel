function j=j_sample_eta(n, jamp, eta)

% j_sample_eta(n, jamp, eta)
% Return a partially symmetric interaction matrix with random Gaussian entries
% and diagonal entries set to zero
% Variance of each element is jamp^2/n, covariance of j_ij j_ji
% eta jamp^2/n

if eta==0
    j = j_sample(n, jamp);
else
    a=sqrt(1-eta^2);
    a=sqrt((1-a)/2);
    b=eta/(2*a);
    j=jamp/sqrt(n)*randn(n);
    j=a*j+b*j';
    j=j-diag(diag(j));
end