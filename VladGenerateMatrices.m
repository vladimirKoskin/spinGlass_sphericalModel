

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
% Initialise
J=zeros(N,N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other 
%%%

% Radius of the N-sphere
R_squared = 1 ;
R = sqrt(R_squared);

% Number of time steps:
N_tot=300000;
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
x(:,1)=x(:,1)*sqrt(N)/sqrt(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Random Matrix
%%(max eigval UNIQUE and positive):


eta=1;

if eta==0
    button = 0;
    while button < 1  
        J=(JAmp/sqrt(N))*randn(N); 
%         J=J-diag(diag(J)) % Non-selfinteracting
        igval=eig(J);
        eigvalmaxPos = find(igval == max(real(eig(J))));
        if imag(igval(eigvalmaxPos)) == 0
            if igval(eigvalmaxPos) < 0.55 & igval(eigvalmaxPos) > 0.0
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
    J=J-diag(diag(J)); % Non-selfinteracting
end;
figure;
plot(real(eig(J)),imag(eig(J)),'o');hold on;grid on;
xlabel('Re(\lambda_J)');
ylabel('Im(\lambda_J)');
hold off;


 %% Generate Random Matrix:
 %%( max eigval COMPLEX CONJUGATE)
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


%% Analytic of the matrix

Jgen = J;
J^(-1)*J;
[V,D] = eig(Jgen);
symJ = (Jgen+Jgen')/2;
skewJ = (Jgen-Jgen')/2;
symJ+skewJ - Jgen;
[Vsym,Dsym] = eig(symJ);
[Vskew,Dskew] = eig(skewJ);
skewJProjectedSym = Vsym'*skewJ*Vsym;
[A,B] = eig(skewJProjectedSym);