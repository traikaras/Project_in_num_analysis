function [B,etri,C]=general_shape(L,hfun)
hfun = +hfun;
%------------------------------------------- setup geometry
    
    node = [                % list of xy "node" coordinates
        0, 1                % outer square
        0, -1
        L, 1
        L, -1 
        L/4, 0.5                % inner square
        L/4, -0.5
        3*L/4, 0.5
        3*L/4, -0.5 ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 4
        3, 4
        1, 3 
        5, 7                % inner square
        7, 8
        6, 8
        5, 6];
%         9, 11                % second inner square
%         11, 12
%         10, 12
%         9, 10] ;
%------------------------------------------- call mesh-gen.
               % uniform "target" edge-lengths
   % B is the coordinate matrix, C the triangle nodes, etri the edge nodes
   [B,etri,C,~] = refine2(node,edge,[],[],hfun) ;
end