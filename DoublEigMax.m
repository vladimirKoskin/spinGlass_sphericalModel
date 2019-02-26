%% Get matrices with max real eigenval double and imaginary part != 0
button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J);
    eigvalmaxPos = find(real(igval) == max(real(eig(J))));
    if imag(igval(eigvalmaxPos)) ~= 0
        button = button + 1;
    end
end

% split in function the following blocks
clear all;

%% Onset Conditions
% dimension:
N=3;
% Radius of the N-sphere
R = sqrt(N);

% correlation paramter:
eta=0.0;
% amplitude
JAmp=1;
% Sparsity of J
Mdensity=0;
% Number of time steps:
N_tot=100000;
% Delta t:
Dt=0.0001;
% Temperature: Increase it!!!
T=0.01;
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
    button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
    igval=eig(J);
    eigvalmaxPos = find(real(igval) == max(real(eig(J))));
    if imag(igval(eigvalmaxPos)) ~= 0
        real(igval(eigvalmaxPos))
        if real(igval(eigvalmaxPos)) < .5 & real(igval(eigvalmaxPos)) > 0.00001
        button = button + 1;
        end
    end
end
 % In case you want sparse:  
 % J=JAmp/sqrt(n)*sprandn(n,n,Mdensity);

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


%% Run multiple times (different IC) for a given Interaction matrix
N_Particles = 5
x=zeros(N*N_Particles,N_tot);
for w=1:N_Particles

    % ICs
var_ICs_X=0.01;
x(:,1) = var_ICs_X*randn(N*N_Particles,1);
M = sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1) = x(:,1)*R/sqrt(M);


% Time steps
var=sqrt(2*Dt*T);

for i=2:N_tot

x(1:3:N_Particles,i)=  x(1:3:N_Particles,i-1)  +  Dt*(-x(1:3:N_Particles,i-1)/N* ( x(1:3:N_Particles,i-1)'*J*x(1:3:N_Particles,i-1)) + J*x(1:3:N_Particles,i-1)) + (eye(N)-1/N*x(1:3:N_Particles,i-1)*x(1:3:N_Particles,i-1)')*randn(N,1)*var;

end

figure(5);hold on
plot3(x(1,:),x(2,:),x(3,:));hold on;grid on
ylim([-(R+0.5) R+0.5]);
xlim([-(R+0.5) R+0.5]);
 title('3-D position'); 
 
end
% Coordinate dynamics
cols=['y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'];
column = 1:N;
figure(98);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,x); hold on
for K = 1 : length(dyncord)
set(dyncord(1),'Color',cols{column(K)});
end
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
plot(x(1,:),x(2,:));
ylim([-(R+1) R+1]);
xlim([-(R+1) R+1]);
hold on
title('2-D position'); hold off


 figure;
plot3(x(1,:),x(2,:),x(3,:));hold on;grid on
ylim([-(R+1) R+1]);
xlim([-(R+1) R+1]);
 title('3-D position'); 


Jgen = J;
J^(-1)*J
[V,D] = eig(Jgen)
symJ = (Jgen+Jgen')/2;
skewJ = (Jgen-Jgen')/2;
symJ+skewJ - Jgen;
[Vsym,Dsym] = eig(symJ)
[Vskew,Dskew] = eig(skewJ)
skewJProjectedSym = Vsym'*skewJ*Vsym;
[A,B] = eig(skewJProjectedSym)
%Endpoint 
x(:, N_tot)/sqrt(x(:, N_tot)'*x(:, N_tot))
J