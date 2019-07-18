close all
clear
clc
%% Setup
n = 2^5; % Number of nodes
L = 1; % Length of the beam
T = 20; % Final time evolution
nt = 200; % Number of time steps
E = 1; % Young's modulus
I = 1; % Area moment of Inertia
mu = 1; % Beam density function

%% Boundary & Initial Conditions
Q_L = -1; % Shear force at pos L
M_L = 0; % Moment at pos L

q = 0; % Load function 
w_0 = 0; % Height of the neutral axis at the first position for all time
w_0_prime = 0; % Derivative of w at the first position for all time

w0 = zeros(2*n+2,1); % Initial position of all nodes at time 0
wp0 = zeros(2*n+2,1); % Initial velocity of all nodes at time 0
wpp0 = zeros(2*n+2,1); % Initial acceleration of all nodes at time 0

B = [1,1,w_0;1,2,w_0_prime;n,3,Q_L;n,4,M_L]; % Boundary condition matrix

%% Form functions
% For the future integration I suppose...
% phi1_bar = @(x) 1-3.*x.^2+2.*x.^3;
% phi2_bar = @(x) x.*(x-1).^2;
% phi3_bar = @(x) 3*x.^2-2*x.^3;
% phi4_bar = @(x) x.^3-x.^2;
% % 
% % Test plot of the functions
% %firstelement = 0:0.01:h;
% %secondelement = h:0.01:h*2;
% %plot(firstelement,phi1_bar((firstelement)./h))
% %hold on 
% %plot(firstelement,phi2_bar((firstelement)./h))
% %plot(firstelement,phi3_bar((firstelement)./h))
% %plot(firstelement,phi4((firstelement)./h))

%% Getting mass and stiffness matrix using it's respective functions
Mass_Matrix = MassMatrix(n,mu);
Stiffness_Matrix = stiffness_matrix(E,I,n);

% 
C = getRestrictionMatrix(B,n);

% Cheat way of getting v_n (Implement the function later)
v_n = zeros(2*n,1);
v_n(end-1) = Q_L;
v_n(end) = M_L;

%% Full extended system
pad_mat = zeros(size(C));
two_by_two = zeros(2,2);
M_l = [Mass_Matrix pad_mat;pad_mat' two_by_two];
S_l = [Stiffness_Matrix C;C' two_by_two];
f = [q+v_n;0;0]; %Extended right hand side

%% Starting from a bent beam
w0 = S_l\f;

%% Time evolution
[W,dt] = time_ev( M_l, S_l, f, w0, wp0 ,wpp0, T, nt );

%% Getting the position values
w = W(1:2:end-2,:);
% and the derivative values
wp = W(2:2:end-2,:);
%% Visualization
show_anim(w,dt,0)

