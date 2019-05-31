% Setup
n = 2^5; % Number of nodes
L = 1; % Length of the beam
tau = 0.2; %time interval
T = 0:tau:10; % Timeline

% Discritization
h = 1/(n-1);
x = 0:h:L; % x position of nodes

% Boundary Conditions
E = 1; % Young's modulus
I = 1; % Area moment of Inertia

mu = 1; % Beam density function
Q_L = -1; % Shear force at pos L
M_L = 0; % Moment at pos L

q = 0; % Load function 

w_0 = 0; % Height of the neutral axis at the first position 
w_0_prime = 0;% Derivative of w at the first position

w_ic0 = zeros(n,1); % initial position of all nodes at time 0
w_prime_ic0 = zeros(n,1); % initial velocity of all nodes at time 0

B = [1,1,w_0;1,2,w_0_prime;n,3,Q_L;n,4,M_L]; % Boundary condition matrix

% Form functions
phi1_bar = @(x) 1-3.*x.^2+2.*x.^3;
phi2_bar = @(x) x.*(x-1).^2;
phi3_bar = @(x) 3*x.^2-2*x.^3;
phi4_bar = @(x) x.^3-x.^2;
% 
% Test plot of the functions
%firstelement = 0:0.01:h;
%secondelement = h:0.01:h*2;
%plot(firstelement,phi1_bar((firstelement)./h))
%hold on 
%plot(firstelement,phi2_bar((firstelement)./h))
%plot(firstelement,phi3_bar((firstelement)./h))
%plot(firstelement,phi4((firstelement)./h))

% Getting mass and stiffness matrix using it's respective functions
Mass_Matrix = MassMatrix(n,mu);
Stiffness_Matrix = stiffness_matrix(E,I,n);

% 
C = getExtendedSystem(B,n);

% Cheat way of getting v_n (Implament the function later)
v_n = zeros(2*n,1);
v_n(end-1) = Q_L;

% Full extended system
pad_mat = zeros(size(C));
two_by_two = zeros(2,2);
M_l = [Mass_Matrix pad_mat;pad_mat' two_by_two];
S_l = [Stiffness_Matrix C;C' two_by_two];
rhs_gamma = [q+v_n;0;0]; 

% Back slash solution of the stationary solution
w_all = S_l\rhs_gamma;

w_clean = w_all(1:2:end-2); 
w_prime = w_all(2:2:end-2);

% figure 
% plot(x,w)
% hold on 
% axis([0 1 -0.8 0.1])

% Newmark Method

gamma = 1/2;
beta = 1/4;

y_Newmark = zeros(2*n+2,1);
y_prime_Newmark = zeros(2*n+2,1);
y_dprime_Newmark = zeros(2*n+2,1);

y_star = zeros(2*n+2,1);
y_prime_star = zeros(2*n+2,1);

y_dprime_Newmark = zeros(2*n+2,1);

for i=2:length(T)
    % Stars
    y_star = y_Newmark + y_prime_Newmark*tau+(0.5-beta)*y_dprime_Newmark*tau^2;
    y_prime_star = y_prime_Newmark+(1-gamma)*y_dprime_Newmark;
    % double prime j+1
    y_dprime_Newmark = (M_l+beta*tau^2*S_l)\(-S_l*y_star);
    % y and y_prime j+1
    y_Newmark = y_star + beta*y_dprime_Newmark*tau^2;
    y_prime_Newmark = y_prime_star + gamma*y_dprime_Newmark*tau;  
    
    figure 
    plot(x,y_Newmark(1:2:end-2))
    hold on 
    axis([0 1 -0.8 0.1])
end

