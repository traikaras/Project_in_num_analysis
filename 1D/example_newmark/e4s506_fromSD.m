%e4s506.m
clear all

t_init = 0; t_final = 3;

t_incr = 0.005;

n = (t_final-t_init)/t_incr;

tt = t_init:t_incr:t_final;

trange = [t_init t_final];

M = [10 0 0;0 20 0;0 0 30];

K = 1e3*[45 -20 -15;-20 45 -25;-15 -25 40];

C = 3e-2*K;

F(1,:) = 0*tt;

F(2,:) = 0*tt;

omega = pi/.29;

F(3,:) = 50*sin(omega*tt);

for j = 1:n
 pulse(j) = sin(omega*tt(j));

 if tt(j) > pi/omega
    pulse(j) = 0;
 end
end

F(3,1:n) = 50*pulse;

x0 = [0 0 0];
xd0= [0 0 0];

[tp,x1] = newmark_fromSD(M,C,K,F,trange,n,x0,xd0);

x1 = 1000*x1;

figure(1), plot(tp,x1(:,1),'k',tp,x1(:,2),'k',tp,x1(:,3),'k',tp,0*tp,'k')
xlabel('Time, s')
ylabel('x, mm')
axis([0 3 -5 10])
hold on
figure(1), plot(tt(1:60),5*pulse(1:60),'.k')
hold off

[t,x2] = ode45(@f10,t_init:t_incr:t_final,[0 0 0 0 0 0]);
x2 = 1000*x2;

d(:,1) = (x1(:,1)-x2(:,1));
d(:,2) = (x1(:,2)-x2(:,2));
d(:,3) = (x1(:,3)-x2(:,3));

figure(2), plot(tp,d(:,1),'k',tp,d(:,2),'k',tp,d(:,3),'k',tp,0*tp,'k')
xlabel('Time, s')
ylabel('Difference, mm')
axis([0 3 -1e-2 1e-2])

%The function f10, required by function ode45, is as follows:


