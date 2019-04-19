function s=sample_noise(n,lengthh,dt,powerspectrum,omegas)

%sample_noise(n,length,dt,powerspectrum,omegas) samples n rows
%of length "lengthh" and with timestep dt
%each row is a random noise process with "powerspectrum" as given over the discrete
%set of angular frequencies "omegas"
%Assume "lengthh" is even here for simplicity, also the lengths
%of the omegas and powerspectrum arrays.

aux=randn(n,lengthh); %sample white noise first
y=fft(aux')';

% go back to discrete power spectrum (divide by dt)
% then interpolate at frequencies required for sample
origlength=length(powerspectrum);
psip=powerspectrum(1:origlength/2+1)/dt;
omip=omegas(1:origlength/2+1);
omq=2*pi/(lengthh*dt)*[0:lengthh/2];
newps=zeros(1,lengthh);
newps(1:lengthh/2+1)=interp1(omip,psip,omq);
newps(end:-1:lengthh/2+2)=newps(2:lengthh/2);
% rescale white noise fourier transform to get required
% correlations
ratio=sqrt(newps/origlength);
yy=y*sparse(1:lengthh,1:lengthh,ratio);
s=ifft(yy')';
