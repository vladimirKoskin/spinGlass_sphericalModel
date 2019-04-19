% example script of how one might solve the dynamical mean field equations
% runs many copies (indicated by n here) of the "one body" (N=1) dynamics,
% initially with zero effective noise chi
% then measures the correlation function (and plots it;
% for the effective noise sampling it actually measures the Fourier transform,
% i.e. the power spectrum) from the second half of the trajectory (to allow the
% one body dynamics to reach stationarity).
% Then samples n copies of effective noise trajectories
% with this power spectrum (times J^2), then repeats.
% Various unchecked assumptions in here:
% - is n large enough and are the trajectories long enough
%   to get reliable statistics on the power spectrum
% - has the power spectrum converged or does one need more iterations
% - is dt small enough (will depend on T)

dt=0.005;
t=200.0;
t0=150;
n=500;
jamp=1.1;
temp=0.1;
x0amp=1.0;
x0=x0amp*randn(n,1);
%x0=x0amp*ones(n,1);
times=[0:dt:t];
effnoise=zeros(n,length(times));
clf;
for i=1:10
    x=one_run_one_body(n,x0,effnoise,dt,temp,t);
    %plot(times,x)
    c=correlation(tanh(x),times,t0,[t0:dt:t],t,dt);
    subplot(2,1,1); hold on; plot([t0:dt:t],c)
    drawnow;
    [ps,om]=power_spectrum(x,2*round(t/(4*dt)),dt);
    pssave(i,:)=ps;
    %size(om)
    effnoise=sample_noise(n,2*round(length(times)/2+0.2),dt,jamp^2*ps,om);
    cc=correlation(effnoise,times,t0,[t0:dt:t],t,dt);
    subplot(2,1,2); hold on; plot([t0:dt:t],cc)
    drawnow;
end