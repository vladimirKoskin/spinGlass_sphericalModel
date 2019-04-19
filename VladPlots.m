PC10 = load('PathCollectionFileDimension10.mat');

% PC10 = struct2cell(PC10);
PC10 = PC10.PathCollection;

path1 = PC10{1,1};
path2 = PC10{1,2};
path3 = PC10{1,3};
figure(7)
for i = 1:50
    atc = atc + autocorr(PC10{1,i}(1,:),10000-1);
end
atc = atc/50
plot(atc)

plot(atc)
plot(autocorr(path1(1,:),10000-1))
plot(autocorr(path1(1,:),10000-1))
plot3(autocorr(path1(1,:)),autocorr(path2(1,:)),autocorr(path3(1,:)))

plot(path1')
% Cross Correlation

xcor = xcorr(  (PC10{1,3}(1,:) -mean(PC10{1,3}(1,:)))/ std(PC10{1,3}(1,:))  ) ;
plot(xcor)

% Normalized cross-correlation
C = normxcorr2(PC10{1,3},PC10{1,3})
figure, surf(C), shading flat

%% MEAN computation:

% <x_i(t)> : can average over temperature but even in absence of
% temperature, can average over ICs
moy = zeros(1,10000)
for j = 1:50
for i = 1:50
    moy = moy + mean(PC10{j,i},1);
    figure(8); hold on
    plot(mean(PC10{j,i},1)) % Nice way to average over the columns (1) or rows (2) or a large matrix, pretty fast!
end
end
hold off
moy = moy/(50*50)

ylabel('Average spin path')
title('Averaged paths over all spins of the system for 50 initial conditions and symmetric random matrices (N=10 , T=0)')
xlabel('Steps (Dt = 0.001)');
plot(moy,'Color','y','LineWidth',3)


%% SPEED computation

for j = 1:50
for i = 1:50
    x = PC10{1,1}
    ParticleSpeed = sqrt(diag((x(:,2:N_tot) - x(:,1:N_tot-1))'*(x(:,2:N_tot) - x(:,1:N_tot-1))));
    ParticleSpeed = x(:,2:N_tot) - x(:,1:N_tot-1);
    moy = moy + mean(PC10{j,i},1);
    figure(8); hold on
    plot(mean(PC10{j,i},1)) % Nice way to average over the columns (1) or rows (2) or a large matrix, pretty fast!
end
end

%% Data Extraction and Analysis

% We have PathCollection{id,im,i} where id is dimensionality index,
% im=matrix index, i=IC index

% We will get aggregates/averages for different dimensionalities separately

% Quantities to compute:
% 1 - Particle speed / distance travelled at each step
% 2 - Particle distance from fixed point(s)
% 3 - Correlation function
% 4 - Average Equilibration Time
% 5 - Spherical constraint holding evolution

Extraction:
FileNum = [2 3 5]

PC2 = load('PathCollectionFileDimension2.mat')
PC2= PC2.PathCollection;
PC3 = load('PathCollectionFileDimension3.mat')
PC3= PC3.PathCollection;
PC5 = load('PathCollectionFileDimension5.mat')
PC5= PC5.PathCollection;
PC10 = load('PathCollectionFileDimension10.mat')
PC10= PC10.PathCollection;

MultiDPathCollection = {PC2 PC3 PC5}
% -- 1 -- %
N_tot = 10000;
NbM = 50;
NbIC = 50;
SpeedCollection = cell(length(FileNum),NbM,NbIC);
ParticleSpeed = zeros(1,N_tot);
for id = 1:length(FileNum)
    
    
for im = 1:50%NbM
for j = 1:50%NbIC
    
   
    x = zeros(FileNum(id),N_tot);
    x(:,:) = cell2mat(MultiDPathCollection{id}(im,j));
% for i = 2:N_tot   
%     ParticleSpeed(i)= sqrt(abs(x(:,i)-x(:,i-1))'*abs(x(:,i)-x(:,i-1)));
% end %%THIS IS NOT OPTIMAL

%-> Look at that

ParticleSpeed = sqrt(diag((x(:,2:N_tot) - x(:,1:N_tot-1))'*(x(:,2:N_tot) - x(:,1:N_tot-1))))
%Yay
SpeedCollection{id,im,j}= ParticleSpeed;
end
end
1+1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for id = 1:length(FileNum)
 
colorstring = 'kbgry';
for im = 1:NbM
for j = 1:NbIC
figure(1);hold on
clr = colorstring(id);

plot(SpeedCollection{id,im,j},'Color',clr);hold on
end
end
end
hold off;hold off;
plot(SpeedCollection{1,1,6})

%%%

% -- 2 -- %

FPDistanceCollection = cell(length(Dimss),NbM,NbIC);
DistFromFP = zeros(1,N_tot);
for id = 1:length(FileNum)
for im = 1:NbM
    J = MatrixCollection{id,im};
    [V,D] = eig(J); %% Get the eigvals and eigvectors of G
    FixedPoint = V(:,find(max(real(eig(J)))))*sqrt(N); % Get the eigenvector associated to the max eigenvalue
for j = i:NbIC
for i = 2:N_tot
    
DistFromFP(i) = sqrt(abs(x(:,i)-stablPoints)'*abs(x(:,i)-stablPoints)); % Compute the distance from eigenvalue at each point
end % i  Iteration
FPDistanceCollection{id,im,j} = DistFromFP;
end % j  Iteration
end % im Iteration
end % id Iteration
%%END2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- 3 -- %
% Tricky, need to think more about that

%%%END3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- 4 -- %

for w=1:NbParticles
equilibrationTime(w) = min(find(distTravelled(w,:) <= 0.00001)) %% This min(find(..)) can be tricky and give errors sometimes
end
figure(9);hold on
plot(equilibrationTime);
mean(equilibrationTime)*Dt
% How to detect equilibration? Either in terms of particle speed vanishing,
% or in terms of distance to fixed point vanishing. Here we get the speed
% stuff but distance might be better (TO COMPARE BOTH RESULTS)

%%%4%%%

% -- 5 -- %

%

C=zeros(N_tot,1);
for v = 1:NbParticles 
    
for i=1:N_tot
C(i)=pathCollection{v}(:,i)'*pathCollection{v}(:,i);
end;
 

figure(41);hold on
plot((0:N_tot-1)*Dt,C);hold on
title('Deviation from sphere'); 

end;
%%%5%%%


%% Plots

col =4; 




%%% 3D Path
for h = 1:NbParticles
    x = pathCollection{h};
    
figure(57);hold on
plot3(x(1,:),x(2,:),x(3,:));hold on
% plot(x(1,:),x(2,:));hold on; grid on
grid on
 ylim([-(R+1.5) R+1.5]);
 xlim([-(R+1.5) R+1.5]);
 title('3-D path (T=0) (Max Eigen Value unique and negative)'); 
 hold off
 
 
 
 
 %%% Coordinate dynamics

figure(98);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,x); hold on
% ylim([-(R+0.5) R+0.5]);
% xlim([0 200]);
% for K = 1 : N
% % set(dyncord(K),'Color',cols(mod(K-1,N)+1));
% set(dyncord(K),'Color',cols(K));
% end;

xlabel('t[unit]');
ylabel('x');
title('Coordinate dynamics'); 
hold off


end
%
 
 


%%%Spherical constraint 

C=zeros(N_tot,1);
for v = 1:NbParticles 
    
for i=1:N_tot
C(i)=pathCollection{v}(:,i)'*pathCollection{v}(:,i);
end;
 

figure(41);hold on
plot((0:N_tot-1)*Dt,C);hold on
title('Deviation from sphere'); 

end;

%%% Equilibration Time
for w=1:NbParticles
equilibrationTime(w) = min(find(distTravelled(w,:) <= 0.00001))
end
figure(9);hold on
plot(equilibrationTime);
mean(equilibrationTime)*Dt
% Tryhard

% stablPointsAlt = V(imag(D) ==0 & real(D)==max(real(eig(J))));
%



% Distance from equilibrium dynamics

figure(91);hold on
% Take into account the transient
dyncord = plot((0:N_tot-1)*Dt,distxeq); hold on
% ylim([-0.5 2*sqrt(N)+0.5]);
% for K = 1 : NbParticles
% % set(dyncord(K),'Color',cols(mod(K-1,N)+1));
% set(dyncord(K),'Color',cols(col));
% end;

%   xlim([990 1000]);
% ylim([5 9]);
% for K = 1 : N
% set(dyncord(K),'Color',cols(mod(K-1,N)+1));
% end
xlabel('t[unit]');
ylabel('x');
title('Relaxation dynamics'); 


%
% Plot Distance traveled at each step for each particle with respect to time

MAV = movmean(distTravelled,3);
figure(102);hold on
% Take into account the transient
col =4;
dyncord = plot((0:N_tot-1)*Dt,distTravelled); hold on
dyncord = plot((0:N_tot-1)*Dt,MAV);hold on
% ylim([-0.5 2*sqrt(N)+0.5]);
% color(col) = 1;
% xlim([0.00001 130]);
% ylim([0.00001 0.005]);
 for K = 1 : NbParticles
% % set(dyncord(K),'Color',cols(mod(K-1,N)+1));
 set(dyncord(K),'Color',cols(col));
end;
xlabel('t[unit]');
ylabel('x');
title('Speed dynamics (distTravelled avg)'); 
% 



figure(66);hold on
for k=1:5
plot(PathCollection{1,k}','Color',colList(k),'LineWidth',2)
end
xlim([0 10000])
title('X-coordinate dynamics (N=3 R=1) Matrix 2');
xlabel('Steps (Dt = 0.001)');
ylabel('Angle (Radian)');