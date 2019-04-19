% Main Code
% Sirio: sirio.belga_fedeli@kcl.ac.uk

% Parallel Computing
% Uncomment the following 2 lines for more machines
%delete(gcp('nocreate'));
%parpool('local');

% split in function the following blocks
clear all;

%% Onset Conditions
% dimension:
N=4;
% Radius of the N-sphere
R = sqrt(N);

% correlation paramter:
eta=0.0;
% amplitude
JAmp=1;
% Sparsity of J
Mdensity=0;
% Number of time steps:
N_tot=1000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.0001;
% Initialise
J=zeros(N,N);
x=zeros(N,N_tot);
% ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*sqrt(N)/sqrt(M);

%% Generate Random Matrix:
if eta==0
    J=JAmp/sqrt(N)*randn(N);   
 % In case you want sparse:  
 % J=JAmp/sqrt(n)*sprandn(n,n,Mdensity);
J=J-diag(diag(J));
else
    a=sqrt(1-eta^2);
    a=sqrt((1-a)/2);
    b=eta/(2*a);    
    J=JAmp/sqrt(N)*randn(N);    
 % In case you want sparse:  
 % J=JAmp/sqrt(n)*sprandn(n,n,Mdensity);
    J=a*J+b*J';
    J=J-diag(diag(J));
end;
figure;
plot(real(eig(J)),imag(eig(J)),'o');hold on;grid on;
xlabel('Re(\lambda_J)');
ylabel('Im(\lambda_J)')
hold off;
%% Time steps

for i=2:N_tot

x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/N* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/N*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);

end;

% alternative Intermediate update: splitting method-

figure;
% Take into account the transient
plot((0:N_tot-1)*Dt,x); hold on;
xlabel('t[unit]');
ylabel('x');
% title('Coordinate dynamics'); hold off



%% Alternative iteration
% Using the lagrangian for future state

%% Evolution of the constraint sum x^2=N
for i=1:N_tot
C(i)=x(:,i)'*x(:,i);
end;

figure;
plot((0:N_tot-1)*Dt,C);hold on
title('Deviation from sphere'); hold off



% 2-D position plot
figure;
plot(x(1,:),x(2,:));hold on
ylim([-(R+1) R+1]);
xlim([-(R+1) R+1]);
hold on
title('2-D position'); 


 figure;
plot3(x(1,:),x(2,:),x(3,:));hold on;grid on
ylim([-(R+1) R+1]);
xlim([-(R+1) R+1]);
 title('3-D position'); hold off



J^(-1)*J
[V,D] = eig(J)

V(1,:)*sqrt(N)
%% Count proba to get H same-real eigenvals
cuntSamReal = 0;
cuntSamImg = 0;
for i=1:10000000
    J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J);
    if real(igval(1))==real(igval(2))
        cuntSamReal=cuntSamReal+1;
    elseif imag(igval(1))==imag(igval(2))
        cuntSamImg = cuntSamImg + 1;
    else
        igval
    end
    
end
eig(J)
imag(igval(1))

%% Get SameReal eigenval matrix
button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J)
    if real(igval(1)) ==  real(igval(2))
        button = button + 1;
    end
end

%% Time steps

for i=2:N_tot

x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/N* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/N*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);

end;

% alternative Intermediate update: splitting method-

figure;
% Take into account the transient
plot((0:N_tot-1)*Dt,x); hold on
xlabel('t[unit]');
ylabel('x');
 title('Coordinate dynamics'); hold off

% 2-D position plot
figure;
plot(x(1,:),x(2,:)); hold on
ylim([-(R+1) R+1]);
xlim([-(R+1) R+1]);
hold on
title('2-D position'); hold off

%% Get OnlyReal eigenval matrix
button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J);
    if imag(eig(J)) ==  0
        button = button + 1;
    end
end

%% The question is: do the largest eigenvalue have a imaginary part of 0 or not? (2 same real but conjugate)

%% Get matrices with max real eigenval unique and imaginary part = 0

button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J);
    eigvalmaxPos = find(igval == max(real(eig(J))));
    if imag(igval(eigvalmaxPos)) == 0
        button = button + 1;
    end
end