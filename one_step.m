function xnew = one_step(n,x,j,dt,temp)

% calculate one timestep dt update from x to xnew at temperature temp
% using specified drift function

noise=sqrt(2*temp*dt)*randn(n,1);
xnew = x + dt*(-x + j*tanh(x)) + noise;