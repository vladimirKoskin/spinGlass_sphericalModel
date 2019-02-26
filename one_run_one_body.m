function x = one_run_one_body(n,x0,effnoise,dt,temp,t)

% one_run_one_body(n,x0,effnoise,dt,temp,t)
% calculate one simulation run over time t, with
% timestep dt at temperature temp, from specified initial condition x0
% and with time history of effective noise as input

steps = round(t/dt);
x=zeros(n,steps+1);
x(:,1)=x0;
for i=1:steps
    x(:,i+1)=one_step_one_body(n,x(:,i),effnoise(:,i),dt,temp);
end