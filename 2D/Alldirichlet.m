close all 
L = 2;
hfun = 0.6;
[B,etri,C] = Rectangle(L,hfun);
  
%% initialization 
rho = 10;
% Lame constants
lambda = 10;
mu = 1;
h = 1; % Thickness of the element
n = length(B); % Number of nodes
Itwo = eye(2); % 2x2 Identity
ntria = length(C); % Number of triangles

%% Boundary Conditions
% Determination of the Dirichlet edges
Moving_nodes = find(B(:,1)==L); % Find nodes on the dir boundary
Nonmoving_nodes = find(B(:,1)==0 | abs(B(:,2))==1);

ind_moving = etri(ismember(etri(:,1),Moving_nodes),:); %Edges with moving nodes
ind_Nonmoving = etri(ismember(etri(:,1),Nonmoving_nodes),:);% Edges with non-moving nodes

% Matrix with nodes of each vertex of the Dir BC
Moving_bound = ind_moving(ismember(ind_moving(:,2),Moving_nodes),:); % Moving edges
Nonmoving_bound = ind_Nonmoving(ismember(ind_Nonmoving(:,2),Nonmoving_nodes),:); % Moving edges
count_boundary = length(Moving_nodes)+length(Nonmoving_nodes); % Amount of Dirichlet BC
% Initial dirichlet BC (values)
u_dir = zeros(2*length(etri),1);
doub = ismember(etri(:,:),Moving_bound);
Moving_row = find(doub(:,1).*doub(:,2));
u_dir(2*Moving_row-1) = 0.1;
%u_dir(7) = 0.2;
count_dir = length(etri);

%C_tilde matrix (C matrix on script)
C_til = zeros(2*n,2*length(etri)); %C from the final system (not FEM C)
% Fill in the matrix with Itwo on specific sites
for i=1:length(etri)
    C_til(2*etri(i,1)-1:2*etri(i,1),2*i-1:2*i) = Itwo;
    C_til(2*etri(i,2)-1:2*etri(i,2),2*i-1:2*i) = Itwo;
end
C_til = sparse(0.5*C_til); 

%% Forces (Tau and gravity)
% Initial external force (gravity)
f = 0*ones(2*ntria,1);
%f(2:2:end-1) = -0.1;

%% Mass and stiffness matrix + Extended system
D_til = 0;
[Me,Se,qe]=extendedsystem(n,B,C,u_dir,count_dir,D_til,C_til,f,rho,lambda,mu,h);

%% Stationary solution
y = Se\qe; 

pStationary = zeros(2*n,2);
pStationary(:,1) = reshape(B',[2*n,1]);
pStationary(:,2) = pStationary(:,1)+y(1:2*n,1);
for i=1:2
    patch('faces',C(:,1:3),'vertices',reshape(pStationary(:,i),[2,n])', ...
        'facecolor','w','edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    title('Stationary solution')
    pause(0.5)
end


% %% Calculating the time evolution
% Time = 20;
% nt = 200;
% dt = Time/nt;
% p = zeros(2*n,nt+1); % Set of all points
% p(:,1) = pStationary(:,1); % Add the initial positions to the matrix
% % Matrix of displacements of all nodes over time
% U = zeros(2*n+2*length(etri),nt+1); 
% Up = zeros(2*n+2*length(etri),1); %Initial velocity
% Upp = zeros(2*n+2*length(etri),1); %Initial acceleration

% %% Solving and plotting for each timestep
% fig = figure;
% grid on
% % x0 = 300;
% % y0 = 300;
% % width = 1500;
% % height = 600;
% % set (gcf, 'position' , [x0, y0, width, height])
% % xlim([-1,11])
% % ylim([-1.5,1.5])
% for i=1:nt
%     %% Solve
%     % Solve using the Newark method
%     [U(:,i+1),Up,Upp] = Newark1step( Me, Se, qe, U(:,i), Up ,Upp,dt);
%     % Store the new positions in a new column of matrix p
%     p(:,i+1) = U(1:2*n,i+1) + p(:,1);
%     
% %     % Stop applying the forces after some time 
% %     if i==1
% %         f =[0;0];
% %         tau_23N=[0;0];
% %     end
% %     % Recalculate the right hand side, for time dependancy
% %     rhs = vertcat(D_til*T + E*f,u_dir);
%     %% Plot
%     if ~ishandle(fig)
%         break
%     end
%     patch('faces',C(:,1:3),'vertices',reshape(p(:,i),[2,n])', ...
%         'facecolor','w', ...
%         'edgecolor',[.2,.2,.2]) ;
%     hold on; axis image off;
%     title(['i=' num2str(i)])
%     pause(0.2)
%     clf;
%     
% end
%     
