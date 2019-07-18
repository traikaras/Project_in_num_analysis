function [B,etri,C]=Rectangle(L,hfun)
hfun = +hfun;
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
               % uniform "target" edge-lengths
   % B is the coordinate matrix, C the triangle nodes, etri the edge nodes
   [B,etri,C,~] = refine2(node,edge,[],[],hfun) ;
end