function plotElementalConstants(u,X,T)

[nelem,nen]=size(T);
p = X'; u = u';
switch nen 
    case 3 %Linear triangles
        t=[T';ones(1,nelem)];
    case 4 %Linear Quadrilaterals (split in two triangles)
        k1=T(:,1); k2=T(:,2); k3=T(:,3); k4=T(:,4); 
        t = [k1 k2 k3;k1 k3 k4]';
        t=[t ; ones(1,2*nelem)];
        u = [u u];
    otherwise
        disp('Error element not implemented')
end
pdesurf(p,t,u) 

