% Try an implementation to the discretization of the model in terms of
% angular coordinates to ensure to stay on the sphere for any Temperature T

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
N_tot=10000;
% Delta t:
Dt=0.001;
% Temperature: Increase it!!!
T=0.000;

N=10;
N_tot=20000
angls = zeros(N-1,N_tot);

%% Transfer ICs

% en.wikipedia.org/wiki/N-sphere    categ "Spherical Coordinates"

id = 2
N=5
 x=zeros(N,N_tot);
 err=zeros(N,50);
for f = 1:50

x(:,1) = ICCollection{id,f};
  x(:,1) = -vF;
for i = 1:N-1
    denomin = sqrt(sum(x(i:N,1).^2)); % = partial R
   angls(i,1) = acos(x(i,1)/denomin);
    
   if i == N-1
       if x(N,1) > 0
           angls(i,1) = acos(x(i,1)/denomin);
       else
%            angls(i,1) =  - acos(x(i,1)/denomin);
           angls(i,1) = 2*pi - acos(x(i,1)/denomin);
       end
   end
   
end
 AngNegVF = angls(:,1)
% OK NOW WE HAVE A FKIN ANGLE VECTOR

% Test reverse process and check if didnt lose info

for i = 1:N
    
    if i < N
    x(i,1) = R*cos(angls(i,1))*prod(sin(angls(1:i-1,1)));
    else
       x(i,1) = R*prod(sin(angls(1:N-1,1))); 
    end
    
end
err(:,f) = x(:,1)- ICCollection{id,f} % Check
end
figure(666);hold on
plot(err) % Error of order E-15 so its alright :p
% title('Difference between initial cartesian coordinate and cartesian coordinates obtained from transforming to spherical coordinates and then back to cartesian for N = 10')
ylabel('Error (order .10^{-16} )');
xlabel('Dimension/component index'); hold off
testx(:,4)
 %Looks not bad

 
 %% test all four quadrants for above
testx = zeros(2,4);
testx(:,1) = [sqrt(2)/2 , sqrt(2)/2]
testx(:,2) = [-sqrt(2)/2 , sqrt(2)/2]
testx(:,3) = [-sqrt(2)/2 , -sqrt(2)/2]
testx(:,4) = [sqrt(2)/2 , -sqrt(2)/2]


%% LETS TRY ITERATION (holy shit)
g = [1 2 3]

NbIC = 10
AngularPath = cell(1,NbIC);

 J = MatrixCollection{1};
 J = MatrixCollection{1,2};
 [U ,I]=eig(J)
 H = zeros(1,N)
 H = 0*H;
 Dt = 0.001;
%
colList = [ 'k' 'b' 'r' 'g' 'm' 'y'];
% FOR 5 IC, 1 matrix in N = 3 :
% FOR 5 IC, 1 matrix in N = 5 :
% FOR 5 IC, 1 matrix in N = 10 : 39.02 seconds
% FOR 5 IC, 1 matrix in N = 20 : 166.25 seconds
% FOR 5 IC, 1 matrix in N = 30 : 
 
tic 
for k = 1:NbIC
x(:,1) = ICCollection{1,k}; % Get initial conditions of dimension id = id and index f
    M=sum(x(:,1).^2);
    % Fix the constraint sum x^2=N for all t
    x(:,1)=x(:,1)*R/sqrt(M);

for i = 1:N-1
    denomin = sqrt(sum(x(i:N,1).^2)); % = partial R
   angls(i,1) = acos(x(i,1)/denomin);
   
   if i == N-1
       if x(N,1) > 0
           angls(i,1) = acos(x(i,1)/denomin);
       else
           angls(i,1) =  2*pi - acos(x(i,1)/denomin);
       end
   end 
end % Get the initial condiions in terms of spherical coordinate (angle vector)

for j = 2:N_tot % Lots of Steps

    % With the idea phi(t+delta) = phi(t) + delta* dot(phi)
   
% for i = 1:N % Loop to recover cartesian version to compute mu nicely with matrix products   
%     
%     if i < N
%     x(i,j-1) = R*cos(angls(i,j-1))*prod(sin(angls(1:i-1,j-1)));
%     else
%        x(i,j-1) = R*prod(sin(angls(:,j-1))); 
%     end   
% end % Output x(t=j) from angls(t=j)
% mu = (x(:,j-1)'*J*x(:,j-1) + H*x(:,j-1) )/R^2;

mu = 0;

for i = 1:N-1
   for u =1:N-1
      
       if i< N
           if u < N
           
            mu = mu + J(i,u)*cos(angls(i,j-1))*cos(angls(u,j-1))*prod(sin(angls(1:i-1,j-1)))*prod(sin(angls(1:u-1,j-1)));
           else 
                 mu = mu + J(i,u)*cos(angls(i,j-1))*prod(sin(angls(1:i-1,j-1)))*prod(sin(angls(1:u-1,j-1)));
           end        
       end
       
      if i ==N 
           
           if u ==N
                 mu = mu + J(N,N)*prod(sin(angls(1:N-1,j-1)))*prod(sin(angls(1:N-1,j-1)));
           else
                 mu = mu + J(i,u)*cos(angls(u,j-1))*prod(sin(angls(1:i-1,j-1)))*prod(sin(angls(1:u-1,j-1)));
           end
           
      end
   end
end
&

%%%
%    for u =1:N-1
%       mu = mu + J(i,u)*cos(angls(i,j-1))*cos(angls(u,j-1))*prod(sin(angls(1:i-1,j-1)))*prod(sin(angls(1:u-1,j-1)));
%    end
% mu = mu + J(i,u)*cos(angls(i,j-1))*prod(sin(angls(1:i-1,j-1)))*prod(sin(angls(1:N-1,j-1)));
%    
% end
% mu = mu + J(N,N)*prod(sin(angls(1:N-1,j-1)))*prod(sin(angls(1:N-1,j-1)));
%%%


phidot = zeros(1,N-1);
for i = 1:N-1 % Loop to update the N-1 angles 
  
    E = 0;
    F = 0;
    for h = 1:(i-1)
        E = E + cot(angls(h,j-1))*phidot(h);      
    end % Output E for this given i and j-1
    
    for h = 1:N-1
        F = F + J(i,h)*cos(angls(h,j-1))*prod(sin(angls(1:h-1,j-1)))  ;  
    end % Output F
    
     F = F + J(i,N)*prod(sin(angls(1:N-1,j-1)));
     F = F / prod(sin(angls(1:i,j-1)));
    
    phidot(i) = cot(angls(i,j-1))*E - F + mu*cot(angls(i,j-1)) + H(i)/(R*prod(sin(angls(1:i,j-1))));
    
    if i == N-1
       angls(i,j) =  mod(angls(i,j-1) + Dt*phidot(i) , 2*pi); 
    else
        angls(i,j) =  mod(angls(i,j-1) + Dt*phidot(i) , pi);
    end
%        angls(i,j) =  angls(i,j-1) + Dt*phidot(i);

end % i iteration (set of new angles computation)

end % j iteration (single path computation)

AngularPath{k} = angls;
figure(8); hold on

plot(angls','Color',colList(k),'LineWidth',2);hold on
%   xlim([0 12000]);
% line([0 10000],[pi/2 pi/2],'LineStyle','--');
%  line([0 10000],[pi/4 pi/4],'LineStyle','--');
%  line([0 10000],[pi*5/4 pi*5/4],'LineStyle','--');
% line([0 10000],[pi*7/4 pi*7/4],'LineStyle','--');
%  line([0 10000],[pi*3/4 pi*3/4],'LineStyle','--');
end % k iteration (multi IC)
toc
title('Angular dynamics (N=3 R=10^5) Matrix 5');
xlabel('Steps (Dt = 0.001)');
ylabel('Angle (Radian)');
hold off; hold off;
figure; hold on
plot(angls');hold on
angls(:,N_tot)



%COmpare with normal iteration

x(:,1) = ICCollection{1,1}
% Time steps
for k = 1:10
    x(:,1) = ICCollection{1,k}
for i=2:N_tot
x(:,i)=  x(:,i-1)  +  Dt*(-x(:,i-1)/(R^2)* ( x(:,i-1)'*J*x(:,i-1)) + J*x(:,i-1)) + (eye(N)-1/(R^2)*x(:,i-1)*x(:,i-1)')*randn(N,1)*sqrt(2*Dt*T);
% distxeq(w,i)= sqrt((x(:,i)-stablPoints)'*(x(:,i)-stablPoints));
% distTravelled(w,i) = sqrt( abs(x(:,i)-x(:,i-1))'*abs(x(:,i)-x(:,i-1)) );
end % End of individual particle path iteration

figure(2); hold on
plot(x(1,:),x(2,:));hold on
xlim([-1 1]);
ylim([-1 1]);
xCart = x;
end

figure(3);hold on
plot(xCart');
plot(xAng');


%recup cartesian coordinates of our angular path
for k = 1:5
x = AngularPath{k};
for i = 1:N
    
    if i < N
    x(i,1) = R*cos(angls(i,1))*prod(sin(angls(1:i-1,1)));
    else
       x(i,1) = R*prod(sin(angls(1:N-1,1))); 
    end
    
end
figure(5);hold on
plot(x');hold on
end
hold off; hold off;
