function PlotNodalField(u,X,T)

[nelem,nen]=size(T);

switch nen
    case 3 % Linear triangles
        p = X';
        t = [T'; ones(1,nelem)];
    case 4 % Linear Quadrilaterals (split in four triangles)
        k1 = T(:,1); k2 = T(:,2); k3 = T(:,3); k4 = T(:,4); k5 = size(X,1) + [1:nelem]';
        XnodesInteriors = (X(k1,:) + X(k2,:) + X(k3,:) + X(k4,:))/4;
        uNodesInteriors = (u(k1) + u(k2) + u(k3) + u(k4))/4;
        u = [u; uNodesInteriors];
        p = [X; XnodesInteriors]';
        t = [k1 k2 k5; k2 k3 k5; k3 k4 k5; k4 k1 k5]';
        t = [t; ones(1,4*nelem)];
    % ... (autres cas pour les différents types d'éléments)
    otherwise
        disp('Error: element not implemented')
        return; % Sortie de la fonction en cas d'erreur
end

pdeplot(p, [], t, 'xydata', u, 'mesh', 'on'); % Utilisation de pdeplot pour afficher la solution et le maillage
