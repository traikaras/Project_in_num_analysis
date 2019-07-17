function [m,s] = localMassStiff(P,rho,lambda,mu,h)
    Itwo = eye(2);
    % Mass matrix
    volB = h*abs(det(P))*0.5;
    m = (kron(ones(3,3),Itwo)+eye(6))*rho*volB/12;
    % Stiffness
    Pinv = P^-1;
    Q = [Pinv(1,1) Pinv(2,1),Pinv(3,1);Pinv(1,2),Pinv(2,2),Pinv(3,2)];
    partE = [1 0 0 0;0 0.5 0.5 0;0 0.5 0.5 0;0 0 0 1];
    E = partE*kron(Q,Itwo);
    lambda_mu = [lambda+2*mu 0 0 lambda; 0 2*mu 0 0; 0 0 2*mu 0;...
        lambda 0 0 lambda+2*mu];
    s = E'*lambda_mu*E*volB;
end