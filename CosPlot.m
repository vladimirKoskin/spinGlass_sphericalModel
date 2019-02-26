a = 1;
t = 0 : 0.01 : 6.3;
A = a*cos(2*t);
B = -a*cos(2*t);

figure(444);
plot(t,B);hold on
plot(t,A);
ylim([-(a+0.1) a+0.1]);