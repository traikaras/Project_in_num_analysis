function Mass=mass_matrix(mu,n)
%n=5;
%mu=1;
h=1/(n-1);
Mass_diag=mu*h/420*sparse(diag([[156,4*h^2], repmat(2*[156,4*h^2], [1, n-2]), [156,4*h^2]]));
Mass_rest=mu*h/420*sparse( diag([22*h, repmat([13*h,0], [1, n-2]), 13*h, -22*h],1)+diag([repmat([54,-3*h^2], [1, n-1])],2)+diag([repmat([-13*h,0], [1, n-2]),-13*h],3));
Mass=Mass_diag+Mass_rest+Mass_rest';
end