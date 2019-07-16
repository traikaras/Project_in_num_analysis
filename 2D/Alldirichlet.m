close all 
L = 10;
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
   hfun = +0.1 ;            % uniform "target" edge-lengths
   [B,etri,C,tnum] = refine2(node,edge,[],[],hfun) ;
   
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
u_dir(2*Moving_row-1) = 0.05;
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
%f(2:2:end) = -0.1;

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
q = E*f;

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
    [m,s] = localMassStiff(palpha,rho);
    
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

%% Stationary solution
y = Se\qe; 

p = zeros(2*n,2);
p(:,1) = reshape(B',[2*n,1]);
p(:,2) = p(:,1)+y(1:2*n,1);
for i=1:2
patch('faces',C(:,1:3),'vertices',reshape(p(:,i),[2,n])', ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    title(['Stationary solution'])
    pause(0.5)
end


%% Calculating the time evolution
Time = 20;
nt = 200;
dt = Time/nt;
p = zeros(2*n,nt+1); % Set of all points
p(:,1) = reshape(B',[2*n,1]); % Add the initial positions to the matrix
% Matrix of displacements of all nodes over time
U = zeros(2*n+2*length(etri),nt+1); 
Up = zeros(2*n+2*length(etri),1); %Initial velocity
Upp = zeros(2*n+2*length(etri),1); %Initial acceleration

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
%     patch('faces',edge(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;
%     title(['i=' num2str(i)])
%     pause(0.2)
%     %clf;
%     
% end
%     
