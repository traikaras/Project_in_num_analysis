function [Stif]= stiffness_matrix(E,I,n)
%n=5;
%E=1;
%I=1;
h=1/(n-1);
i=1:2*n;
Stif_diag=(E*I/h^3)*sparse(diag([[12,4*h^2], repmat(2*[12,4*h^2], [1, n-2]), [12,4*h^2]]));
Stif_rest=(E*I/h^3)*sparse( diag([6*h, repmat([-6*h,0], [1, n-2]), -6*h, -6*h],1)+diag([repmat([-12,2*h^2], [1, n-1])],2)+diag([repmat([6*h,0], [1, n-2]),6*h],3));
Stif=Stif_diag+Stif_rest+Stif_rest';
end