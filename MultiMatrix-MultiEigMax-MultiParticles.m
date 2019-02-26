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


%%

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
N_tot=200000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.0;
% Initialise
J=zeros(N,N);
x=zeros(N,N_tot);
% ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*sqrt(N)/sqrt(M);

%% Generate Random Matrix (max eigval unique and positive):

maxEigenvalWanted = 0.2; % will be between value and value-0.1


delt = 0.1;

matrixPerEigMax = 5;

MatrixCollection = cell(1,110);
maxEigenvalWanted = 0.2;
cunty=1;
for maxEigenvalWanted = 0.2:1:1.2
    
    for f=1:matrixPerEigMax
        
        
if eta==0
    button = 0;
while button < 1  
 J=JAmp/sqrt(N)*randn(N); 
 
    igval=eig(J);
    eigvalmaxPos = find(igval == max(real(eig(J))));
    if imag(igval(eigvalmaxPos)) == 0
        if igval < maxEigenvalWanted & igval > maxEigenvalWanted-delt
        button = button + 1;
        end
    end
end  


end
  MatrixCollection{cunty} = J;
  cunty = cunty+1;
    end
end
 % In case you want sparse:  
 % J=JAmp/sqrt(n)*sprandn(n,n,Mdensity);

% else
%     a=sqrt(1-eta^2);
%     a=sqrt((1-a)/2);
%     b=eta/(2*a);    
%     J=JAmp/sqrt(N)*randn(N);    
%  % In case you want sparse:  
%  % J=JAmp/sqrt(n)*sprandn(n,n,Mdensity);
%     J=a*J+b*J';
%     J=J-diag(diag(J));
% end;

for i =1:110
  plot(real(eig(MatrixCollection{i})),imag(eig(MatrixCollection{i})),'o');hold on 
end


figure;
plot(real(eig(J)),imag(eig(J)),'o');hold on;grid on;
xlabel('Re(\lambda_J)');
ylabel('Im(\lambda_J)')
hold off;


%% Run multiple times (different IC) for a given Interaction matrix
nbParticles = 2;
nbMatrix = 10;
PathCollection = cell(nbMatrix,nbParticles);


cols=['y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'];
column = 1:N;

for mat=1:nbMatrix

    J = MatrixCollection{mat};
for k=1:nbParticles
    
    
   for i=2:N_tot

    x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/N* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/N*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);

    end;
    
    PathCollection{mat,k} =  x;
    
    
    % Plot coordinates dynamics
   figure(98);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,x); hold on
ylim([-(R+0.5) R+0.5]);
xlim([0 80]);
set(dyncord(1:3),'Color',cols(fix(mat/matrixPerEigMax)+1));
xlabel('t[unit]');
ylabel('x');
title('Coordinate dynamics'); 
legend;

end

end

plot((0:N_tot-1)*Dt,PathCollection{1,1});




cols=['y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'];
column = 1:N;
var=sqrt(2*Dt*T);
for w=1:100
    
    % ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*R/sqrt(M);


% Time steps

for i=2:N_tot

x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/N* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/N*x(:,i-1)*x(:,i-1)')*randn(N,1)*var;

end

figure(56);hold on
plot3(x(1,:),x(2,:),x(3,:));hold on;grid off
color(1) = 0.5;
ylim([-(R+1.5) R+1.5]);
xlim([-(R+1.5) R+1.5]);
 title('3-D path (T=0) (Max Eigen Value unique and negative)é'); 
 
% Coordinate dynamics

figure(98);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,x); hold on
ylim([-(R+0.5) R+0.5]);
xlim([180 200]);
for K = 1 : N
set(dyncord(K),'Color',cols(mod(K-1,N)+1));
end
xlabel('t[unit]');
ylabel('x');
title('Coordinate dynamics'); 
end
