
%% Onset Conditions

% Radius of the N-sphere
R_squared = 1 ;
R = sqrt(R_squared);

% Number of time steps:
N_tot=10000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.000;
% Initialise
x=zeros(N,N_tot);

% ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*sqrt(N)/sqrt(M);


%% Run multiple times (different IC) for a given Interaction matrix

NbParticles=10;
cols=['y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'];
column = 1:N;
var=sqrt(2*Dt*T);
distxeq = zeros(NbParticles,N_tot);
distTravelled = zeros(NbParticles,N_tot);
pathCollection = cell(NbParticles,1);

%    J = symJ;
%     J=skewJ;
%    

Transmat = J-eye(N)*(x(:,1)'*J*x(:,1)/N);
[V,D] = eig(Transmat);
stablPoints = V(:,find(max(real(eig(J)))))*sqrt(N);

for w=1:NbParticles

    % ICs
var_ICs_X=0.01;
x(:,1)=var_ICs_X*randn(N,1);
M=sum(x(:,1).^2);
% Fix the constraint sum x^2=N for all t
x(:,1)=x(:,1)*R/sqrt(M);



i=2
w=1
% Time steps

for i=2:N_tot
x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/(R^2)* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/(R^2)*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);
distxeq(w,i)= sqrt((x(:,i)-stablPoints)'*(x(:,i)-stablPoints));
distTravelled(w,i) = sqrt( abs(x(:,i)-x(:,i-1))'*abs(x(:,i)-x(:,i-1)) );


end;

pathCollection{w} = x;

% equilibrationTime(w) = min(find(distTravelled(w,:) < 0.000001));

% for i=1:N_tot
% PrD(i)=expRotAxis'*x(:,i);
% end;
% figure(39);hold on
% plot((0:N_tot-1)*Dt,PrD);hold on
% title('Projection of trajectory on skew-symmetric eigenvector related to eigval= 0'); hold off
% hold off
 
end;
distTravelled(:,1)=distTravelled(:,2);