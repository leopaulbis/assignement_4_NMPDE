function plotElementalConstants1D(rho,X,T)

hold on
nOfElements = size(T,1);
for e=1:nOfElements
    Te=T(e,[1,2]);
    Xe=X(Te);
    plot(Xe,rho(e)*[1,1],'r-')
end
hold off