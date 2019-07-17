function [ edge ] = getEdge( coord, B, list_edges )
%GETEDGE gets a point in 2D space and return the indices of the edge that is
%the closest to the point, i.e., the row in the matrix of edges. The output
%is the array of indices of the boundary condition vector that should be
%modified by the correspondant force or displacement, depending on if its
%Dirichlet of Neumann
%   The function recieves as parameter the 1x2 coordinate coord, the matrix
%   of coordinates of all nodes B, and the matrix containing the nodes of
%   each edge in its columns list edges.
%   The coordinate should be in the boundary

un = unique(list_edges(:));
b = B(un,:);
diff = sqrt(sum((coord-b).^2,2)); % Find distances between nodes
[~,idx] = sort(diff); % Sort distances
c = [b(idx(1),:); b(idx(2),:)]; % Choose two closest nodes
ind_nodes = find(ismember(B,c,'rows')==1);
edge_node = find(ismember(list_edges,ind_nodes','rows')==1);
edge = [2*edge_node-1,2*edge_node];

end

