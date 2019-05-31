function [tp,x] = newmark_fromSD(M,C,K,F,tspan,n,x0,xd0)

% M, C, K are matrices multiplying xddot, xdot and x repectively

% F is column vector of exciations. x0, xd0 are initial x0 and xd vectors

% tspan = [t_initial t_final]; n is tspan/t_increment

dt = (tspan(2)-tspan(1))/n;

tp(1) = tspan(1);

x(:,1) = x0';

xd(:,1) = xd0';

gamma = 1/2; beta = 1/4;

A = (1/(beta*dt^2))*M+(gamma/(beta*dt))*C+K; 
invA = inv(A);

xdd(:,1) = inv(M)*(F(:,1)-C*xd(:,1)-K*x(:,1));

for i = 1:n

 B = (F(:,i+1)+ M*((1/(beta*dt^2))*x(:,i)+(1/(beta*dt))*xd(:,i)+(1/(2*beta)-1)*xdd(:,i))+C*((gamma/(beta*dt))*x(:,i)+(gamma/beta-1)*xd(:,i)+(gamma/beta-2)*(dt/2)*xdd(:,i)));

 x(:,i+1) = invA*B;

 xdd(:,i+1) = (1/(beta*dt^2))*(x(:,i+1)-x(:,i))-(1/(beta*dt))*xd(:,i)-((1/(2*beta))-1)*xdd(:,i);

 xd(:,i+1) = xd(:,i)+(1-gamma)*dt*xdd(:,i)+gamma*dt*xdd(:,i+1);

 tp(i+1) = tp(i)+dt;

end

x = x';