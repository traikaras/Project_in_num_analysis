close all 
L = 2;
node = [                % list of xy "node" coordinates
        0, 1                % outer square
        0, -1
        L, 1
        L, -1 
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 3                % outer square 
        3, 4
        2, 4
        2, 1 
        ] ;
%------------------------------------------- call mesh-gen.
   hfun = +1.5 ;            % uniform "target" edge-lengths
   [B,etri,C,tnum] = refine2(node,edge,[],[],hfun) ;
   
%% initialization 
rho = 10;
% Lame constants
lambda = 1;
mu = 1;
h = 1; % Thickness of the element
n = length(B); % Number of nodes
Itwo = eye(2); % 2x2 Identity
ntria = length(C); % Number of triangles
%% Boundary Conditions
% Determination of the Dirichlet edges
dir_nodes = find(B(:,1)==0); % Find nodes on the dir boundary
ind_dir = etri(ismember(etri(:,1),dir_nodes),:); % Edges with dir nodes
% Matrix with nodes of each vertex of the Dir BC
dir_bound = ind_dir(ismember(ind_dir(:,2),dir_nodes),:); % Dir edges
count_dir = length(dir_bound); % Amount of Dirichlet BC
% Initial dirichlet BC (values)
u_dir = zeros(length(dir_bound(:)),1);

%C_tilde matrix (C matrix on script)
C_til = zeros(2*n,2*count_dir); %C from the final system (not FEM C)
% Fill in the matrix with Itwo on specific sites
for i=1:count_dir
    C_til(2*dir_bound(i,1)-1:2*dir_bound(i,1),2*i-1:2*i) = Itwo;
    C_til(2*dir_bound(i,2)-1:2*dir_bound(i,2),2*i-1:2*i) = Itwo;
end
C_til = sparse(0.5*C_til); 

% Determination of Neumann edges
neu_nodes = find(B(:,1)~=0 | abs(B(:,2))==1);
ind_neu = etri(ismember(etri(:,1),neu_nodes),:);
% Matrix with nodes of each vertex of Neu BC
neu_bound = ind_neu(ismember(ind_neu(:,2),neu_nodes),:);
count_neu = length(neu_bound); %Amount of Neumann BC

%D_tilde matrix (Multiplied with known taus)
D_til = zeros(2*n,2*count_neu);% Equivalent of C_til for the neumann BC
% Fill in the matrix with Itwo on specific sites
for i=1:count_neu
    D_til(2*neu_bound(i,1)-1:2*neu_bound(i,1),2*i-1:2*i) = Itwo;
    D_til(2*neu_bound(i,2)-1:2*neu_bound(i,2),2*i-1:2*i) = Itwo;
end
D_til = sparse(0.5*D_til);
%% Forces (Tau and gravity)
% Initial Neumann BC
T = zeros(length(D_til(1,:)),1);
%T(42) = -1;
% Initial external force (gravity)
f = -0*ones(2*ntria,1);
%f(2:2:end) = -0.1;

% Apply a tau on the top right corner (comment later)
% ind_corner = find(neu_bound==find(B(:,1)== L & B(:,2)==1));
% bla = neu_bound([3 4],:);
% T(2*ind_corner(B(bla(:,2),2)==1))=-2000;
%% RHS
% E matrix: Matrix multiplied by f. Defines all the triangles
E = zeros(2*n,2*ntria);
% Fill in the matrix with Itwo on specific sites
for i=1:ntria
    E(2*C(i,1)-1:2*C(i,1) , 2*i-1:2*i) = Itwo;
    E(2*C(i,2)-1:2*C(i,2) , 2*i-1:2*i) = Itwo;
    E(2*C(i,3)-1:2*C(i,3) , 2*i-1:2*i) = Itwo;
end
E = sparse(E/3);
% Right Hand Side
q = D_til*T + E*f;
%% Global Stiffness and Mass Matrix
% We form the diagonal and upper part separately. Summ upper and transpose
% in the final result taking advantage of symmetry
MG = zeros(2*n,2*n);
MDiag = zeros(2*n,2*n);
SG = zeros(2*n,2*n);
SDiag = zeros(2*n,2*n);

for alpha=1:ntria
    % Calculate the local Mass and stiff matrix for each triangle
    palpha = B(C(alpha,:),:)';
    palpha = vertcat(palpha,ones(1,3));
    [m,s] = localMassStiff(palpha,rho,lambda,mu,h);
    
    %MASS
    % Positions in Global matrix of local values
    MGpos1 = 2*C(alpha,1)-1;
    MGpos2 = 2*C(alpha,2)-1;
    MGpos3 = 2*C(alpha,3)-1;
    % Fill in the diagonal
    MDiag(MGpos1:MGpos1+1,MGpos1:MGpos1+1)=MDiag(MGpos1:MGpos1+1,MGpos1:MGpos1+1)+m(1:2,1:2);
    MDiag(MGpos2:MGpos2+1,MGpos2:MGpos2+1)=MDiag(MGpos2:MGpos2+1,MGpos2:MGpos2+1)+m(3:4,3:4);
    MDiag(MGpos3:MGpos3+1,MGpos3:MGpos3+1)=MDiag(MGpos3:MGpos3+1,MGpos3:MGpos3+1)+m(5:6,5:6);
    % Fill in the upper part
    MG(MGpos1:MGpos1+1,MGpos2:MGpos2+1)=MG(MGpos1:MGpos1+1,MGpos2:MGpos2+1)+m(1:2,3:4);
    MG(MGpos1:MGpos1+1,MGpos3:MGpos3+1)=MG(MGpos1:MGpos1+1,MGpos3:MGpos3+1)+m(1:2,5:6);
    MG(MGpos2:MGpos2+1,MGpos3:MGpos3+1)=MG(MGpos2:MGpos2+1,MGpos3:MGpos3+1)+m(3:4,5:6);
    %STIFFNESS
    % Positions in Global matrix of local values
    SGpos1 = 2*C(alpha,1)-1;
    SGpos2 = 2*C(alpha,2)-1;
    SGpos3 = 2*C(alpha,3)-1;
    % Fill in the diagonal
    SDiag(SGpos1:SGpos1+1,SGpos1:SGpos1+1)=SDiag(SGpos1:SGpos1+1,SGpos1:SGpos1+1)+s(1:2,1:2);
    SDiag(SGpos2:SGpos2+1,SGpos2:SGpos2+1)=SDiag(SGpos2:SGpos2+1,SGpos2:SGpos2+1)+s(3:4,3:4);
    SDiag(SGpos3:SGpos3+1,SGpos3:SGpos3+1)=SDiag(SGpos3:SGpos3+1,SGpos3:SGpos3+1)+s(5:6,5:6);
    % Fill in the upper part
    SG(SGpos1:SGpos1+1,SGpos2:SGpos2+1)=SG(SGpos1:SGpos1+1,SGpos2:SGpos2+1)+s(1:2,3:4);
    SG(SGpos1:SGpos1+1,SGpos3:SGpos3+1)=SG(SGpos1:SGpos1+1,SGpos3:SGpos3+1)+s(1:2,5:6);
    SG(SGpos2:SGpos2+1,SGpos3:SGpos3+1)=SG(SGpos2:SGpos2+1,SGpos3:SGpos3+1)+s(3:4,5:6);
end
% Sum the diagonal, upper and lower parts
SG = SG'+SG+SDiag;
MG = MG'+MG+MDiag;

%% Build extended system
Zbig = zeros(size(C_til)); %Zero block of size C_til
Zsmall = zeros(2*count_dir,2*count_dir);%Lower right corner zero block
% MASS
Me = sparse([MG Zbig;Zbig' Zsmall]);
% STIFFNESS
Se = sparse([SG C_til; C_til' Zsmall]);
% RHS
qe = [q;u_dir];

%% Calculating the time evolution
Time = 5;
nt = 50;
dt = Time/nt;
p = zeros(2*n,nt+1); % Set of all points
p(:,1) = reshape(B',[2*n,1]); % Add the initial positions to the matrix
% Matrix of displacements of all nodes over time
U = zeros(2*n+2*count_dir,nt+1); 
Up = zeros(2*n+2*count_dir,1); %Initial velocity
Upp = zeros(2*n+2*count_dir,1); %Initial acceleration

%% Solving and plotting for each timestep
fig = figure;
grid on
% x0 = 300;
% y0 = 300;
% width = 1500;
% height = 600;
% set (gcf, 'position' , [x0, y0, width, height])
% xlim([-1,11])
% ylim([-1.5,1.5])
for i=1:nt
    %% Solve
    % Solve using the Newark method
    [U(:,i+1),Up,Upp] = Newark1step( Me, Se, qe, U(:,i), Up ,Upp,dt);
    % Store the new positions in a new column of matrix p
    p(:,i+1) = U(1:2*n,i+1) + p(:,1);
    
%     % Stop applying the forces after some time 
%     if i==1
%         f =[0;0];
%         tau_23N=[0;0];
%     end
%     % Recalculate the right hand side, for time dependancy
%     rhs = vertcat(D_til*T + E*f,u_dir);
    %% Plot
    if ~ishandle(fig)
        break
    end
    patch('faces',C(:,1:3),'vertices',reshape(p(:,i),[2,n])', ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['i=' num2str(i)])
    pause(dt)
    %clf;
    
end
    