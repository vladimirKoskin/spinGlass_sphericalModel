function c=correlation(x,times,t1,t2,t,dt)

% correlation(x,times,t1,t2,t,dt)
% calculates correlation between x-values at times t1 and t2
% where t1 and t2 can be arrays
% and t, dt are total time and time step of data points in x
%
% average is also done over different runs in x (if more than one)

ind1=1+round(t1/dt);
ind2=1+round(t2/dt);
    
% case of a single run
if length(size(x))==2
    n=size(x);
    n=n(1); % number of "spins"
    c=x(:,ind1)'*x(:,ind2)/n;
% case of multiple runs
else
    n=size(x);
    runs=n(1);
    n=n(2);
    c=0;
    for i=1:runs
        c=c+reshape(x(i,:,ind1),length(ind1),n)*...
            reshape(x(i,:,ind2),n,length(ind2));
    end
    c=c/(n*runs);
end
    
    