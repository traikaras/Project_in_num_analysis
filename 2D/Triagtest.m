    node = [                % list of xy "node" coordinates
        0, 1                % outer square
        0, -1
        10, 1
        10, -1 
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 3                % outer square 
        3, 4
        2, 4
        2, 1 
        ] ;

%------------------------------------------- call mesh-gen.
    %hfun = +1.7 ;            % uniform "target" edge-lengths
    hfun = +0.7 ;
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun) ;

%------------------------------------------- draw tria-mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    temp = B(etri);
    for i=1:length(temp)
        text(B(i,1),B(i,2),num2str(i))
    end