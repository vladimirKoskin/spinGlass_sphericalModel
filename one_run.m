function x = one_run(n,x0,j,dt,temp,t)

% calculate one simulation run over time t, with
% timestep dt at temperature temp, from specified initial condition x0
% and with given interaction matrix j

steps = round(t/dt);
x=zeros(n,steps+1);
x(:,1)=x0;
for i=1:steps
    x(:,i+1)=one_step(n,x(:,i),j,dt,temp);
end