% example script of how to do a simulation run and plot x_i(tau)

dt=0.001;
t=2000.0;
n=20;
eta=0.0;
jamp=2/(1+eta);
temp=0.4;
x0amp=1.0;
%j=j_sample_eta(n,jamp,eta);
%x0=x0amp*randn(n,1);
x=one_run(n,x0,j,dt,temp,t);
times=[0:dt:t];
plot(times,x)