function xnew = one_step_one_body(n,x,effnoise,dt,temp)

% calculate one timestep dt update from x to xnew at temperature temp
% using specified drift function and a vector of effective noise
% variables

noise=sqrt(2*temp*dt)*randn(n,1);
xnew = x + dt*(-x+effnoise) + noise;