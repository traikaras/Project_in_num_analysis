function [Me,Se,qe]=extendedsystem(n,B,C,u_dir,count_dir,D_til,C_til,f,rho,lambda,mu,h)
ntria = length(C);
Itwo = eye(2); % 2x2 Identity

% RHS
% E matrix: Matrix multiplied by f. Defines all the triangles
E = zeros(2*n,2*ntria);
% Fill in the matrix with Itwo on specific sites
for i=1:ntria
    E(2*C(i,1)-1:2*C(i,1) , 2*i-1:2*i) = Itwo;
    E(2*C(i,2)-1:2*C(i,2) , 2*i-1:2*i) = Itwo;
    E(2*C(i,3)-1:2*C(i,3) , 2*i-1:2*i) = Itwo;
end
E = sparse(E/3);

% Global Stiffness and Mass Matrix
% We form the diagonal and upper part separately. Summ upper and transpose
% in the final result taking advantage of symmetry
MG = zeros(2*n,2*n);
MDiag = zeros(2*n,2*n);
SG = zeros(2*n,2*n);
SDiag = zeros(2*n,2*n);

for alpha=1:ntria
    % Calculate the local Mass and Stiff matrix for each triangle
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

% Build extended system
Zbig = zeros(size(C_til)); %Zero block of size C_til
Zsmall = zeros(2*count_dir,2*count_dir);%Lower right corner zero block
% MASS
Me = sparse([MG Zbig;Zbig' Zsmall]);
% STIFFNESS
Se = sparse([SG C_til; C_til' Zsmall]);
% RHS
T = zeros(length(D_til(1,:)),1);
q = D_til*T + E*f;
qe = [q;u_dir];
end
