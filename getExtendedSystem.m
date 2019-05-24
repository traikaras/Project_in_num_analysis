function [ C ] = getExtendedSystem(B,M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = length(M(:,1))/2; % We might want to just take in n directly 

% B = [1,1,0;1,2,0;n,3,0;n,4,1];
f1 = find(B(:,2)==1); % What if there are more than one with index 1 ?
f2 = find(B(:,2)==2);
e = full(ind2vec((2*B(f1,1)-1)'));
delta1 = zeros(2*n-length(e(:,1)),length(e(1,:)));
e = vertcat(e,delta1);
d = full(ind2vec((2*B(f2,1))'));
delta2 = zeros(2*n-length(d(:,1)),length(d(1,:)));
d = vertcat(d,delta2); % Compines the vectors into a matrix [d;delta2]
%d = wextend('addrow','zpd',d,2*n-length(d(:,1)),'d');
C = [e d];

end

