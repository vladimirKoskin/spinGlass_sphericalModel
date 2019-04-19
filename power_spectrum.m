function [p,freqs]=power_spectrum(x, length, dt)

%power_spectrum(x, length) takes the last "length" columns of x
%and computes their fft, and from this the power spectrum
%(averaged over rows)
%Multiplication by dt ensures that the power spectrum is normalized
%as it would be in continuous time
%(i.e. so that for small dt the Fourier sums become integrals)

y=fft(x(:,end-length+1:end)')';
z=abs(y).^2;
p=mean(z)*dt;
freqs=2*pi/(length*dt)*[0:length-1];