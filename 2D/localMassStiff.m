function [m,s] = localMassStiff(P,rho,lambda,mu,h)
    Itwo = eye(2);
    % Mass matrix
    volB = h*abs(det(P))*0.5;
    m = (kron(ones(3,3),Itwo)+eye(6))*rho*volB/12;
    % Stiffness
    Q = P^-1;
    partE = [1 0 0 0;0 0.5 0.5 0;0 0.5 0.5 0;0 0 0 1];
    E = partE*kron(Q(1:2,:),Itwo);
    lambda_mu = [lambda+2*mu 0 0 lambda; 0 2*mu 0 0; 0 0 2*mu 0; lambda 0 0 lambda+2*mu];
    s = volB*E'*lambda_mu*E;
end