function [ C ] = getRestrictionMatrix( B,n)
%GETEXTENDEDSYSTEM Computes the restriction matrix C of the system given
%the matrix of initial conditions B. The output is a 2*n+2 times number of
%B.C. matrix with ones on the position corresponding to nodes with B.C. 
%   B is the matrix of initial conditions. Its first column represent the
%   number of the node where there is a boundary. The second column the
%   type of boundary (1,2,3 or 4) and the last column the numerical value.
%   n is the number of nodes

f1 = find(B(:,2)==1); % Find B.C. of type 1 (Dirichlet)
f2 = find(B(:,2)==2); % Find B.C. of type 2 (Neumann)

%% Put a 1 in the position of the node which has BC of type 1
e = full(ind2vec((2*B(f1,1)-1)'));
delta1 = zeros(2*n-length(e(:,1)),length(e(1,:)));
e = vertcat(e,delta1);
%% Put a 1 in the position of the node which has BC of type 2
d = full(ind2vec((2*B(f2,1))'));
delta2 = zeros(2*n-length(d(:,1)),length(d(1,:)));
d = vertcat(d,delta2);
%% Concatenate the corresponding matrices
C = [e d];
end

