%% Lame values
rho = 10;
lambda = 1;
mu = 1;
f = [0;0];
h = 1;
p_1 = [0;3];
p_2 = [1;0];
p_3 = [-1;0];

%% Boundary conditions
u_1D = 0.5*(p_1+p_2);
u_3D = 0.5*(p_1+p_3);
tau_23N =[0;0];
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
C = 0.5*[Itwo Itwo;Itwo ztwo;ztwo Itwo];
S_extend = zeros(10);
S_extend(1:6,:)=[S C];
S_extend(7:10,1:6)=C';
% Right hand side
rhs = vertcat(T_23*tau_23N+V*f,u_1D,u_3D);

%% Calculating the time evolution
u0 = zeros(10,1);
u0(1:6,1) = [p_1;p_2;p_3];
wp0 = zeros(10,1);
wpp0 = zeros(10,1);
[U,u_prime] = time_ev( M_extend, S_extend, rhs, u0, wp0 ,wpp0, 20, 200 );

%% Plotting up the solution
u1 = U(1:2,:);
u2 = U(3:4,:);
u3 = U(5:6,:);
X = [u1(1,:);u2(1,:);u3(1,:)];
Y = [u1(2,:);u2(2,:);u3(2,:)];
fig = figure;
axis([-2 2 -2 5])
for i=1:100
    if ~ishandle(fig)
        break
    end
    patch(X(:,i),Y(:,i),'r')

    pause(0.5)
    %clf;
end


%show_anim(U(:end-4,:),50);