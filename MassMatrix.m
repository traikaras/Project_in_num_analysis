function [Mass] = MassMatrix(n,mu)
%Assebly of mass matrix for the wave equation
%   n is the number of nodes
%   mu is the density of the beam (assumed to be constant
h = 1/(n-1);
Mass_diagonal = mu*h/420*sparse(diag([156,4*h^2, repmat(2*[156,4*h^2], [1, n-2]...
    ), 156,4*h^2]));

Mass_offdiagonals = mu*h/420*sparse(diag([22*h, repmat([13*h,0], [1, n-2])...
    , 13*h, -22*h],1)+diag([repmat([54,-3*h^2], [1, n-1])],2)+...
    diag([repmat([-13*h,0], [1, n-2]),-13*h],3));

Mass = Mass_diagonal+Mass_offdiagonals+Mass_offdiagonals';
end

