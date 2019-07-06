    
node = [                % list of xy "node" coordinates
        0, 1                % outer square
        0, -1
        2, 1
        2, -1 
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 3                % outer square 
        3, 4
        2, 4
        2, 1 
        ] ;
%------------------------------------------- call mesh-gen.
   hfun = +1.7 ;            % uniform "target" edge-lengths
   [B,etri,C,tnum] = refine2(node,edge,[],[],hfun) ;
   
%% initialization 
rho = 10;
lambda = 1;
mu = 1;
f = [0;0];
h = 1;
n = length(B);
% BC

%% 
MG = zeros(2*n,2*n);
MDiag = zeros(2*n,2*n);

for alpha=1:length(C)
    palpha = B(C(alpha,:),:)';
    palpha = vertcat(palpha,ones(1,3));
    [m,s] = localMassStiff(palpha,rho);
    %for j=1:3
%     MDiag(C(alpha,1),1) = MDiag{C(alpha,1),1}+ m(1:2,1:2);
%     MDiag(C(alpha,2),1) = MDiag(C(alpha,2),1)+ m(3:4,3:4);
%     MDiag(C(alpha,3),1) = MDiag(C(alpha,3),1)+ m(5:6,5:6);
%     %end
%     MG(C(alpha,1),C(alpha,2)) = MG(C(alpha,1),1)+ m(1:2,3:4);
%     MG(C(alpha,1),C(alpha,3)) = MG(C(alpha,1),C(alpha,3))+ m(1:2,5:6);
%     MG(C(alpha,2),C(alpha,2)) = MG(C(alpha,2),C(alpha,2))+ m(3:4,5:6);
    for i=1:3 % We only need to loop over node ones element
        two = 1;
        three = 1;
        
        MGpos = 2*C(alpha,i)-1;
        % if the node is not on the boundary the part beloning only to it self should go into the global matrix
        MDiag(MGpos:MGpos+1,MGpos:MGpos+1)=MDiag(MGpos:MGpos+1,MGpos:MGpos+1)+m(2*i-1:2*i,2*i-1:2*i);

        if two % this is true if the first node beloging to this element was true
            Mtwo = 2*C(alpha,2)-1;
            MG(Mtwo:Mtwo+1,MGpos:MGpos+1)=MG(Mtwo:Mtwo+1,MGpos:MGpos+1)+m(2*i-1:2*i,3:4);
        end

        if three % this is true if the second node beloging to this element was true
            Mthree = 2*C(alpha,3)-1;
            MG(Mthree:Mthree+1,MGpos:MGpos+1)=MG(Mthree:Mthree+1,MGpos:MGpos+1)+m(2*i-1:2*i,5:6);
        end

        % setting the boolean values for the first and second node
        if i==1
            two=0;
        elseif i==2
            three=0;
        end
    end
end
MG = MG'+MG+MDiag;



