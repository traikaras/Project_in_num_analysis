close all
clear
clc
%% Setup
n = 2^5; % Number of nodes
L = 1; % Length of the beam
T = 20; % Final time evolution
nt = 200; % Number of time steps
tstop = T; %Time when force is set to zero
E = 1; % Young's modulus
I = 1; % Area moment of Inertia
mu = 1; % Beam density function
begin_st_state = 0; % Boolean to determine if we start from steady state
mov = 0; % Set to one to save video in avi format
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
if begin_st_state
    w0 = S_l\f;
end

%% Time evolution
[W,dt] = time_ev( M_l, S_l, f, w0, wp0 ,wpp0, T, nt, tstop );

%% Getting the position values
w = W(1:2:end-2,:);
% and the derivative values
wp = W(2:2:end-2,:);
%% Visualization
show_anim(w,dt,mov)

