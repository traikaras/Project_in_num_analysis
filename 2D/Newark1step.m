function [W,Wp,Wpp] = Newark1step( M, S, f, w, wp ,wpp,dt)
%TIME_EV takes the system of differential equations described with mass
%matrix M, stifness matrix S and right hand side f, and evolves it one
%timestep
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

A = M + beta*dt^2*S; % Matrix for solving the linear system en each step

% MOVE INTO A Function
%% Step a) Calculation of intermediate values
% Equation 7
w_s = w + wp*dt +(0.5-beta)*wpp*dt^2;
% Equation 8
wp_s = wp + (1-gamma)*wpp*dt;
%% Step b) Solving the linear system for second derivative
p = f-S*w_s;
Wpp = A\p;
%% Step c) Getting w and w prime with Equation 9
W = w_s + beta*dt^2*Wpp;
Wp = wp_s + gamma*dt*Wpp;
end