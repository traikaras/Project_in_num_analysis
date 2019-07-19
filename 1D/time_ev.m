function [ W, dt ] = time_ev( M, S, f, w0, wp0 ,wpp0, T, nt,t_stop )
%TIME_EV takes the system of differential equations described with mass
%matrix M, stifness matrix S and right hand side f, and evolves it from
%time zero up to time T, with nt timesteps. 
%   This scripts implements Newmark's method for time evolution of
%   differential equations for the case of a beam under the action of an
%   external force. All the dynamics of the system is given in matrices M
%   and S. We assume that the system has no dissipative terms (D=0) and
%   that matrices M,S and f are time independent.
%   w0 is a vector of positions of each node of the beam at time 0 as a
%   column vector.
%   wp0 is a vector of first derivative of positions of each node of the
%   beam at time 0 as a column vector
%   wpp0 is a vector of second derivative of positions of each node of the
%   beam at time 0 as a column vector

% Parameters for the weighted mean of the derivatives needed for the method
beta = 0.25;
gamma = 0.5;
dt = T/nt; % size of timestep (constant)

% Initialize the output with fixed dimension
n = length(w0);
W = zeros(n,nt+1);
W(:,1) = w0; % The first state is the initial state

% The derivatives in previous step initialized as the initial values given.
wp = wp0; 
wpp = wpp0;
A = M + beta*dt^2*S; % Matrix for solving the linear system en each step

for i=2:nt+1
    % MOVE INTO A Function
    %% Step a) Calculation of intermediate values
    % Equation 7
    w_s = W(:,i-1) + wp*dt +(0.5-beta)*wpp*dt^2;
    % Equation 8
    wp_s = wp + (1-gamma)*wpp*dt;
    %% Step b) Solving the linear system for second derivative
    p = f-S*w_s;
    wpp = A\p;
    %% Step c) Getting w and w prime with Equation 9
    W(:,i) = w_s + beta*dt^2*wpp;
    wp = wp_s + gamma*dt*wpp;
    
    %-------------------------------TEST--------------------------------
    % Time dependent Q_L 
    % if Q_L starts at zero this slowly increases a downward force 
%     if i*dt<5
%         f(end-3)=f(end-3)-i*0.05;
%     end

    % At time point t_stop it sets Q_L as 0
    if i*dt==t_stop
        f(end-3)=0;
    end
end

