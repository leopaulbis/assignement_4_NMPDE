function plotDerivative1D(X,T,u)
grado = size(T,2)-1;

xi = linspace(-1,1,11)';

hold on
if grado == 1
    unos=ones(11,1);
    N=[(1-xi)/2 (1+xi)/2];
    dNdxi=[unos*(-1/2) unos*(1/2)];
elseif grado == 2
    N = [(xi-1).*xi/2   (1+xi).*xi/2    (1-xi).*(1+xi)];
    dNdxi=[xi-1/2 xi+1/2 -2*xi];
end
x_interp = []; du_interp = [];
for ielem = 1:size(T,1)
    Te = T(ielem,:);
    Xe = X(Te,:);
    ue = u(Te);
    J=dNdxi*Xe; %Jacobian at Gauss points
    dNdx = diag(1./J)*dNdxi; %derivatives with respect to x
    x_interp = N*Xe;
    du_interp = dNdx*ue;
    plot(x_interp,du_interp,'k');
end
    
hold off