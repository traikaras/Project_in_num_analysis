%% Lame values
rho = 1;
lambda = 10^4;
mu = 10^4;
f = [100;50];
h = 2;
p_1 = [0;0];
p_2 = [0;1];
p_3 = [0.5;0.5];

%% Boundary conditions
u_1D = p_1;
u_3D = p_3;
tau_23N =[1;1];
%% Setup
Itwo = eye(2);
ztwo = zeros(2);
P = vertcat([p_1 p_2 p_3],ones(1,3));
Q = P^-1;
stupidE = [1 0 0 0;0 0.5 0.5 0;0 0.5 0.5 0;0 0 0 0.5];
E = stupidE*kron(Q(1:2,:),Itwo);
volB = h*abs(det(P))/2;
% T's
l_12 = norm(p_1-p_2);
l_23 = norm(p_2-p_3);
l_31 = norm(p_1-p_3);
T_12 = l_12*h/2*[Itwo;Itwo;ztwo];
T_23 = l_23*h/2*[ztwo;Itwo;Itwo];
T_31 = l_31*h/2*[Itwo;ztwo;Itwo];

V = volB/3*[Itwo;Itwo;Itwo];
%% Mass and stiffness matrix 
M = ones(6)+eye(6)*rho*volB/12;
lambda_mu = [lambda+2*mu 0 0 lambda; 0 2*mu 0 0; 0 0 2*mu 0; lambda 0 0 lambda+2*mu];
S = E'*lambda_mu*E;

%% Extended system
% Mass matrix
M_extend = zeros(10);
M_extend(1:6,1:6)=M;
% Stiffness matrix
S_extend = zeros(10);
S_extend(1:6,1:10)=[S -T_12 -T_31];
S_extend(7:8,1:2)=Itwo;
S_extend(9:10,5:6)=Itwo;
% Right hand side
rhs = vertcat(T_23*tau_23N+V*f,u_1D,u_3D);

w0 = zeros(6,1);
wp0=zeros(6,1);
wpp0=zeros(6,1);
[U,u_prime] = time_ev( M_extend, S_extend, rhs, w0, wp0 ,wpp0, 20, 50 );

show_anim(U,50);