
[I,M] = eig(J)
vF = I(:,1)

lastpos = angls(:,N_tot)

lastposX 
for k=1:10
    
path = AngularPath{k};
diff1 = path(:,N_tot) - AngNegVF
diff2 = path(:,N_tot) - AngPosVF
err = min(diff1'*diff1 ,diff2'*diff2)
plot(err);hold on
end
path = AngularPath{9};
path(:,N_tot) - AngPosVF
path(:,N_tot) -AngNegVF
AngPosVF


%% Angular -> Cartesian
lastposX = zeros(1,N)
for i = 1:N
    
    if i < N
    lastposX(i) = R*cos(angls(i,N_tot))*prod(sin(angls(1:i-1,N_tot)));
    else
       lastposX(i) = R*prod(sin(angls(1:N-1,N_tot))); 
    end
    
end

%% Cartesian -> angular

PathCollection{5,1}(:,N_tot)
lastposAng = zeros(1,N-1);
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
