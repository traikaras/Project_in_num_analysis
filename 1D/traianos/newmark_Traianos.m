function [tp,w] = newmark_Traianos(M_l,K_l,S_l,rhs_gamma,tau,dt,w_0,w_0_prime)
% M_l,K_l,S_l are matrices multiplying xddot, xdot and x repectively
% F is column vector of exciations. x0, xd0 are initial x0 and xd vectors
% tspan = [t_initial t_final]; n is tspan/t_increment
tp(1) = 0;
w(:,1) = w_0;
wd(:,1) = w_0_prime;
gamma = 1/2; beta = 1/4;
A = (1/(beta*dt^2))*M_l+(gamma/(beta*dt))*K_l+S_l; 
% for general case wdd(:,1) = M_l\(rhs_gamma(:,1)-K_l*wd(:,1)-S_l*w(:,1));
wdd(:,1)=zeros(size(wd));

for i = 1:tau
 B = (rhs_gamma(:)+...
 M_l*((1/(beta*dt^2))*w(:,i)+(1/(beta*dt))*wd(:,i)+ ...
(1/(2*beta)-1)*wdd(:,i))+K_l*((gamma/(beta*dt))*w(:,i)+...
(gamma/beta-1)*wd(:,i)+(gamma/beta-2)*(dt/2)*wdd(:,i)));

 w(:,i+1) = A\B;
 
 wdd(:,i+1) = (1/(beta*dt^2))*(w(:,i+1)-w(:,i))...
-(1/(beta*dt))*wd(:,i)-((1/(2*beta))-1)*wdd(:,i);

 wd(:,i+1) = wd(:,i)+(1-gamma)*dt*wdd(:,i)+gamma*dt*wdd(:,i+1);
 tp(i+1) = tp(i)+dt;

end
w = w';