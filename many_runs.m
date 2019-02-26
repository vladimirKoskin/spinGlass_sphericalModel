function x = many_runs(n,x0,j,dt,temp,t,nruns)

% usage: x = many_runs(n,x0,j,dt,temp,t,nruns)
% calculate nruns simulation run over time t, with
% timestep dt at temperature temp, from specified initial condition x0
% and with given interaction matrix j

steps = round(t/dt);
x=zeros(nruns,n,steps+1);
for run=1:nruns
    x(run,:,1)=x0;
    for i=1:steps
        x(run,:,i+1)=one_step(n,x(run,:,i)',j,dt,temp);
    end
end