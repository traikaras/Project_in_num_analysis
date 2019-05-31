function yprime = f10(t,y)

m = [10 0 0;0 20 0;0 0 30];

k = 1e3*[45 -20 -15;-20 45 -25;-15 -25 40];

c = 3e-2*k; 

f =[0; 0; 50];

omega = pi/.29;

pulse = sin(omega*t);

if t > pi/omega
 pulse = 0;
end

A = [zeros(3,3) eye(3,3); -m\k -m\c];
b = [zeros(3,1); m\f*pulse];

yprime = A*y+b;
end

