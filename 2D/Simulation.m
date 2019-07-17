% This script solves the motion of a beam with FEM. The beam has length L
% and thickness 2. Is centered in y=0 and the right extremum is at x=0. The
% boundary conditions are either Dirichlet (we set by hand the
% displacement) or Neumann (we set by hand the external force. In each
% simulation the Dirichlet and Neumann borders are chosen.The beam has an
% internal weight force represented by f and it acts over each element. 
% The triangulization id donde automatically with the mesh2d package.

close all 
%% initialization 
rho = 1000; % Density
% Lame constants
lambda = 1000;
mu = 100000;
h = 1; % Thickness of the element
L = 50; %Length of Beam
% 'magic number' for the triangulation. The lower, the more triangles.
hfun = +0.3; 
Itwo = eye(2); % 2x2 Identity
Time = 50; % Final time for time evolution
nt = 500; % Number of time steps
dt = Time/nt; % Size of timestep
delay = dt/10;
begin_st_state = 1; % Boolean to determine if we start from steady state

%% Triangulization
[B,etri,C] = Rectangle(L,hfun);
% Constants from triangulization
n = length(B); % Number of nodes
ntria = length(C); % Number of triangles
%% Forces
% Initial external force (gravity)
f = zeros(2*ntria,1);
f(2:2:end) = -0.015; % Y component of the weight
coord = [L,1]; % Specific place
tau = [0,-1]; % Force aplied
%% Boundary Conditions (Should be a function)
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
% Apply a tau on a specific place
edge = getEdge(coord,B,neu_bound); % Indices of T to be modified
T(edge) = tau; %Replace
%% Build Extended system
% Extended mass, extended stiffness and extended right hand side, and E
% respectively
[Me,Se,qe,E] = extendedsystem(n,B,C,u_dir,count_dir,D_til,C_til,f,rho,lambda,mu,h,T);
%% Calculating the time evolution
p = zeros(2*n,nt+1); % Set of all points
p(:,1) = reshape(B',[2*n,1]); % Add the initial positions to the matrix
% Matrix of displacements of all nodes over time
U = zeros(2*n+2*count_dir,nt+1); 
Up = zeros(2*n+2*count_dir,1); %Initial velocity
Upp = zeros(2*n+2*count_dir,1); %Initial acceleration

%% Starting from stationary solution
if begin_st_state
    u_stat = Se\qe;
    U(:,1) = u_stat;
    pstart = p(:,1) + u_stat(1:2*n);
end
%% Solving and plotting for each timestep


for i=1:nt
    %% Solve
    % Solve using the Newark method
    [U(:,i+1),Up,Upp] = Newark1step( Me, Se, qe, U(:,i), Up ,Upp,dt);
    % Store the new positions in a new column of matrix p
    p(:,i+1) = U(1:2*n,i+1) + p(:,1);
    
    % Stop applying the forces after some time 
    if i==20
        T = zeros(length(D_til(1,:)),1);
        f = zeros(2*ntria,1);
        f(2:2:end) = 0.05;
    end
    % Recalculate the right hand side, for time dependancy
    qe = vertcat(D_til*T + E*f,u_dir);
 
end

%% Plot
p(:,1) = pstart;
fig = figure;
grid on
for j=1:nt
    
    if ~ishandle(fig)
        break
    end
    clf;
    patch('faces',C(:,1:3),'vertices',reshape(p(:,j),[2,n])', ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
%     patch('faces',edge(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;
    title(['Time=' num2str(j*dt)])
    pause(delay)
    
end
    