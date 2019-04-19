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

% split in function the following blocks
clear all;

%% Onset Conditions
% dimension:
N=15;
% Radius of the N-sphere
R_squared = 1 ;
R = sqrt(R_squared);

% correlation paramter:
eta=0.0;
% amplitude
JAmp=1;
% Sparsity of J
Mdensity=0;
% Number of time steps:
N_tot=200000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.000;
% Initialise
J=zeros(N,N);
x=zeros(N,N_tot);

% ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*R/sqrt(M);

%% Generate Random Matrix (max eigval unique and positive):


if eta==0
    button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
 
    igval=eig(J);
    eigvalmaxPos = find(igval == max(real(eig(J))));
    if imag(igval(eigvalmaxPos)) == 0
        if igval(eigvalmaxPos) < 0.55 & igval(eigvalmaxPos) > 0.5
        button = button + 1;
        end;
    end;
end;  
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
ylabel('Im(\lambda_J)');
hold off;


%% Run multiple times (different IC) for a given Interaction matrix

NbParticles=10;
cols=['y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'];
column = 1:N;
var=sqrt(2*Dt*T);
distxeq = zeros(NbParticles,N_tot);
distTravelled = zeros(NbParticles,N_tot);


for w=1:NbParticles
    
    % ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*R/sqrt(M);


Transmat = J-eye(N)*(x(:,1)'*J*x(:,1)/N);
[V,D] = eig(Transmat);
stablPoints = V(:,find(max(real(eig(J)))))*sqrt(N);

% Time steps
i=1
for i=2:N_tot

x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/R_squared* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/N*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);
distxeq(w,i)= sqrt((x(:,i)-stablPoints)'*(x(:,i)-stablPoints));
distTravelled(w,i) = sqrt(abs(x(:,i)-x(:,i-1))'*abs(x(:,i)-x(:,i-1)));

end;

end;

%Decouple Nb of dimensions and radius of sphere
%% Plots

col =4;

% Tryhard

% stablPointsAlt = V(imag(D) ==0 & real(D)==max(real(eig(J))));
%
% figure(56);hold on
%plot3(x(1,:),x(2,:),x(3,:));hold on;grid off
% color(col) = 1;
% ylim([-(R+1.5) R+1.5]);
% xlim([-(R+1.5) R+1.5]);
%  title('3-D path (T=0) (Max Eigen Value unique and negative)'); 
 
% Coordinate dynamics

figure(98);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,x); hold on
% ylim([-(R+0.5) R+0.5]);
% xlim([0 200]);
for K = 1 : N
% set(dyncord(K),'Color',cols(mod(K-1,N)+1));
set(dyncord(K),'Color',cols(col));
end;
color(col) = 1;
xlabel('t[unit]');
ylabel('x');
title('Coordinate dynamics'); 

% Distance from equilibrium dynamics
figure(99);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,distxeq); hold on
% ylim([-0.5 2*sqrt(N)+0.5]);
for K = 1 : NbParticles
% set(dyncord(K),'Color',cols(mod(K-1,N)+1));
set(dyncord(K),'Color',cols(col));
end;

% xlim([1 40]);
% for K = 1 : N
% set(dyncord(K),'Color',cols(mod(K-1,N)+1));
% end
xlabel('t[unit]');
ylabel('x');
title('Relaxation dynamics'); 


%
% Plot Distance traveled at each step for each particle with respect to time

figure(101);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,distTravelled); hold on
% ylim([-0.5 2*sqrt(N)+0.5]);
color(col) = 1;
%  xlim([1 10]);
for K = 1 : NbParticles
% set(dyncord(K),'Color',cols(mod(K-1,N)+1));
set(dyncord(K),'Color',cols(col));
end;
xlabel('t[unit]');
ylabel('x');
title('Speed dynamics (distTravelled)'); 




% !!!! Implement a list of matrices with A = cell(1,100); to record all
% paths x_i and to store multiple matrices

%% Alternative colors

figure(2);hold on
% Take into account the transient
plot((0:N_tot-1)*Dt,x,'Color',[x,0,1]); hold on;
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

end

J^(-1)*J
[V,D] = eig(J)
