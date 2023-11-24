function theReferenceElement = createReferenceElement(degree,typeOfElement)
%theReferenceElement = createReferenceElement(degree,typeOfElement)
%Creates a struct with the information of the reference element, for
%triangle discretizations

switch degree
    case 1
        if typeOfElement==1 %TRI
            %element pour faire la quadrature (calcul de l'integrale dans
            %l'élément)
            z = [0.5 0; 0 0.5; 0.5 0.5];%integration points
            w = [1 1 1]/6;%weight

            xi = z(:,1); eta=z(:,2);
            N = [1-xi-eta,   xi,   eta];%fonction nodale dans l'élément de référence (N1,N2,N3) évaluée sur les points de gauss car on n'as besoin de les connaitres que aux points de gauss 
            Nxi  = [-ones(size(xi))   ones(size(xi))    zeros(size(xi))]; %derivées par rapport à xi de N1,N2 et N3 évaluées aux points de gauss (on prend size(xi) car on veut un tableau avec les coordonnées. Ici, les dérivées étant constantes (valent 1 ou zeros) d'où le ones ou zero)
            Neta = [-ones(size(xi))   zeros(size(xi))   ones(size(xi))];%dérivéees par rapport à eta de N1,N2 et N3 évaluées aux points de gauss 

            nodesCoord = [0,0;1,0;0,1];%coordonées des noeuds dans l'élément dans lequel les points d'intégrations doivent de trouvées
        else
            error('Element not implemented yet')
        end
    case 2
        
         [z,w]=getTriangleQuadrature(5);
         xi = z(:,1); eta=z(:,2);

         N1=zeros(numel(xi),1);%contain the values of the nodal function N1 
         % évaluated in all the integration points
         N2=zeros(numel(xi),1);
         N3=zeros(numel(xi),1);
         N5=zeros(numel(xi),1);
         N4=zeros(numel(xi),1);
         N6=zeros(numel(xi),1);
         
         N1xi=zeros(numel(xi),1);
         N2xi=zeros(numel(xi),1);
         N3xi=zeros(numel(xi),1);
         N4xi=zeros(numel(xi),1);
         N5xi=zeros(numel(xi),1);
         N6xi=zeros(numel(xi),1);
            
         N1eta=zeros(numel(xi),1);
         N2eta=zeros(numel(xi),1);
         N3eta=zeros(numel(xi),1);
         N4eta=zeros(numel(xi),1);
         N5eta=zeros(numel(xi),1);
         N6eta=zeros(numel(xi),1);

         for i=1:numel(xi)
            N1(i) = (1 - xi(i) - eta(i)) * (1 - 2 * xi(i) -2*eta(i));
            N2(i) = 2 * xi(i)^2- xi(i);
            N3(i) = 2 * eta(i)^2 - eta(i);
            N4(i) = 4 * xi(i) * (1 - xi(i) - eta(i));
            N5(i) = 4 * xi(i) * eta(i);
            N6(i) = 4 * eta(i) * (1 - xi(i) - eta(i));

            N1xi(i)=-3+4*xi(i)+4*eta(i);
            N2xi(i)=4*xi(i)-1;
            N4xi(i)=4-8*xi(i)-4*eta(i);
            N5xi(i)=4*eta(i);
            N6xi(i)=-4*eta(i);

            N1eta(i)=4*xi(i)+4*eta(i)-3;
            N3eta(i)=4*eta(i)-1;
            N4eta(i)=-4*xi(i);
            N5eta(i)=4*xi(i);
            N6eta(i)=4-4*xi(i)-8*eta(i);
         end

         N=horzcat(N1,N2,N3,N4,N5,N6); %nodal functions
         Nxi = [N1xi,N2xi,N3xi,N4xi,N5xi,N6xi];
         Neta = [N1eta,N2eta,N3eta,N4eta,N5eta,N6eta];

        nodesCoord = [0 0;1 0;0 1;1/2 0;1/2 1/2; 0 1/2];
        
end
if typeOfElement==1
    typeElement='TRI';
else
    typeElement='QUA';
end
theReferenceElement=struct('IPweights',w,'IPcoord',z,'N',N,'Nxi',Nxi,'Neta',Neta,'type',typeElement,'nodesCoord',nodesCoord);
