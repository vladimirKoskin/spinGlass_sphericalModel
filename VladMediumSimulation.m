%% Parameters

clear all;

%% Onset Conditions

% Radius of the N-sphere
R = 1;
% amplitude
JAmp=1;
% Sparsity of J
Mdensity=0;
% Number of time steps:
N_tot=50000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.000;


%% MODAFUKIN RUN DA SIMULATION
NbM = 5; % Nb of matrices per collection
Dimss = [100];
PathCollection = cell(NbM,NbIC); % Crossing fingers it won't end up being 50 Go data

% PathCollection = zeros(NbM,NbIC,Dimss(id),N_tot);

PathCollection = cell(NbM,NbIC);
 N=Dimss(id);
 id=2;
for im = 1:NbM % 'im' denotes the matrix index (Iteration over each matrix for a given dimensionality)
% im=1
    J=MatrixCollection{im};
for w= 1:NbIC % Iteration over each IC for a given matrix in a given dimensionality
% w=1
    x=zeros(N,N_tot);
x(:,1) = ICCollection{w}; % We have same ICs for different matrices of the same dimensionality, for easiness and better 
%--comparativity
    M=sum(x(:,1).^2);
    % Fix the constraint sum x^2=N for all t
    x(:,1)=x(:,1)*R/sqrt(M);
% Time steps

for i=2:N_tot
x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/(R^2)* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/(R^2)*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);
% distxeq(w,i)= sqrt((x(:,i)-stablPoints)'*(x(:,i)-stablPoints));
% distTravelled(w,i) = sqrt( abs(x(:,i)-x(:,i-1))'*abs(x(:,i)-x(:,i-1)) );
end % End of individual particle path iteration

% distTravelled(:,1)=distTravelled(:,2);
PathCollection{im,w} = x;

% equilibrationTime(w) = min(find(distTravelled(w,:) < 0.000001));
 
end %END w iteration (over ICs for a given dimensionality)
end %END im iteration (over matrices for a given dimensionality)
% NamePC = 'PathCollectionDimesion'+string(Dimss(id))+'TEST.mat';
% save(char(NamePC),'PathCollection','-v7.3');
% clear PathCollection;

figure(5);hold on
plot(PathCollection{1,4}')
xlim([0 5000])
title('Coordinate dynamics (N=3 R=1) Matrix 1')
xlabel('Steps (dt = 0.001)');
ylabel('Xi');
%Check:
PathCollection{1,1,1}
PathCollection{length(Dimss),NbM,NbIC}
PC = PathCollection;