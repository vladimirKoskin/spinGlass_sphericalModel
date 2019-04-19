% example script of how to do multiple simulation runs;
% a correlation function is computed at the end

dt=0.01;
t=500.0;
n=50;
n_runs=10;
jamp=1.5;
temp=0.1;
x0amp=1.0;
j=j_sample_eta(n,jamp,0);
%x0=x0amp*ones(n,1);
x0=x0amp*randn(n,1);
x=many_runs(n,x0,j,dt,temp,t,n_runs);
times=[0:dt:t];
t1=t/2;
t2=[t1:dt:t];
c=correlation(x,times,t1,t2,t,dt);
plot(t2-t1,c)