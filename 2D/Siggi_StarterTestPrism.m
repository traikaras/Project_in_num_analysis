%% Lame values
rho = 10;
lambda = 1;
mu = 1;
f = [0;0];
h = 1;
p1 = [0;3];
p2 = [1;0];
p3 = [-1;0];

%% Boundary conditions
u_1D = 0.5*(p1+p2);
u_3D = 0.5*(p1+p3);
tau_23N =[0;-1];
%% Setup
Itwo = eye(2);
ztwo = zeros(2);
P = vertcat([p1 p2 p3],ones(1,3));
Q = P^-1;
stupidE = [1 0 0 0;0 0.5 0.5 0;0 0.5 0.5 0;0 0 0 0.5];
E = stupidE*kron(Q(1:2,:),Itwo);
volB = h*abs(det(P))/2;
% T's
l_12 = norm(p1-p2);
l_23 = norm(p2-p3);
l_31 = norm(p1-p3);
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
C = 0.5*[Itwo Itwo;Itwo ztwo;ztwo Itwo];
S_extend = zeros(10);
S_extend(1:6,:)=[S C];
S_extend(7:10,1:6)=C';
% Right hand side
rhs = vertcat(T_23*tau_23N+V*f,u_1D,u_3D);

%% Calculating the time evolution
T = 20;
nt = 200;
dt = T/nt;
p = zeros(6,nt+1);
p(:,1) = vertcat(p1,p2,p3);
U = zeros(10,nt+1);
Up = zeros(10,1);
Upp = zeros(10,1);
for i=1:nt
    [U(:,i+1),Up,Upp] = Newark1step( M_extend, S_extend, rhs, U(:,i), Up ,Upp,dt);
    p1 = U(1:2,i+1)+p1;
    p2 = U(3:4,i+1)+p2;
    p3 = U(5:6,i+1)+p3;
    l_12 = norm(p1-p2);
    l_23 = norm(p2-p3);
    l_31 = norm(p1-p3);
    T_12 = l_12*h/2*[Itwo;Itwo;ztwo];
    T_23 = l_23*h/2*[ztwo;Itwo;Itwo];
    T_31 = l_31*h/2*[Itwo;ztwo;Itwo];
    p(:,i+1) = vertcat(p1,p2,p3);
    P = vertcat([p1 p2 p3],ones(1,3));
    volB = h*abs(det(P))/2;
    Q = P^-1;
    E = stupidE*kron(Q(1:2,:),Itwo);
    V = volB/3*[Itwo;Itwo;Itwo];
    M = ones(6)+eye(6)*rho*volB/12;
    M_extend(1:6,1:6)=M;
    S = E'*lambda_mu*E;
    S_extend(1:6,:)=[S C];
    S_extend(7:10,1:6)=C';
    rhs = vertcat(T_23*tau_23N+V*f,u_1D,u_3D);
end


%% Plotting up the solution
% u1 = p(1:2,:);
% u2 = p(3:4,:);
% u3 = p(5:6,:);
% X = [u1(1,:);u2(1,:);u3(1,:)];
% Y = [u1(2,:);u2(2,:);u3(2,:)];
Px = p(1:2:5,:);
Py = p(2:2:6,:);
fig = figure;
%axis([-2 2 -2 5])
for i=1:100
    if ~ishandle(fig)
        break
    end
    patch(Px(:,i),Py(:,i),'r')

    pause(0.5)
    %clf;
end


%show_anim(U(:end-4,:),50);