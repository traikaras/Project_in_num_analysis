% Setup
n = 2^3; % Number of nodes
L = 1; % Length of the beam
T = 0:10; % Time ?

% Discritization
tau = 1./T; 
h = 1/(n-1);
x = 0:h:L; % x position of nodes

% Boundary Conditions
E = 1; % Young's modulus
I = 1; % Area moment of Inertia

mu = 1; % Beam density function
Q_L = 1; % Shear force at pos L
M_L = 1; % Moment at pos L

q = 0; % Load function 

w_0 = 0; % Height of the neutral axis at the first position 
w_0_prime = 0;% Derivative of w at the first position

B = [1,1,w_0;1,2,w_0_prime;n,3,Q_L;n,4,M_L]; % Boundary condition matrix

% Form functions
phi1_bar = @(x) 1-3.*x.^2+2.*x.^3;
phi2_bar = @(x) x.*(x-1).^2;
phi3_bar = @(x) 3*x.^2-2*x.^3;
phi4_bar = @(x) x.^3-x.^2;
% 
% % Test plot of the functions
% firstelement = 0:0.01:h;
% secondelement = h:0.01:h*2;
% plot(firstelement,phi1_bar((firstelement)./h))
% hold on 
% plot(firstelement,phi2_bar((firstelement)./h))
% plot(firstelement,phi3_bar((firstelement)./h))
% plot(firstelement,phi4((firstelement)./h))

% Getting mass and stiffness matrix using it's respective functions
Mass_Matrix = MassMatrix(n,mu);
Stiffness_Matrix = stiffness_matrix(E,I,n);

% 
C = getExtendedSystem(B,n);
% Cheat way of getting v_n (Implament the function later)
v_n = zeros(2*n,1);
v_n(end-1) = 1;

% Full extended system
pad_mat = zeros(size(C));
two_by_two = zeros(2,2);
M_l = [Mass_Matrix pad_mat;pad_mat' two_by_two];
S_l = [Stiffness_Matrix C;C' two_by_two];
rhs_gamma = [q+v_n;0;0]; 

% Newark Method


