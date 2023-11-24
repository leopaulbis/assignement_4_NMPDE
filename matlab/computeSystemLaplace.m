function [K,f]=computeSystemLaplace(X,T,theReferenceElement,sourceTermFunction)

IPweights = theReferenceElement.IPweights;
IPcoord = theReferenceElement.IPcoord;
N=theReferenceElement.N;%function evaluated on the gauss point
Nxi=theReferenceElement.Nxi; %dérivative with respect to xi on the gauss point
Neta=theReferenceElement.Neta;%same with respect to eta

%number of nodes =number of row of the matrix which contains the coordinate of the nodes
nOfNodes = size(X,1); 
nOfElements =size(T,1);%number of element=nb of row in the connectivity matrix
%initialize a sparse matrix of size n0nodes*n0nodes whith room for up to 9*nonodes 
K=spalloc(2*nOfNodes,2*nOfNodes,18*nOfNodes);
f=zeros(2*nOfNodes,1);%initialise the vector of the second member
n1=nOfNodes;
K1=K(1:n1,1:n1);
K2=K(n1+1:end,n1+1:end);
K3=K(1:n1,n1+1:end);
K4=K(n1+1:end,1:n1);

f2=f(n1+1:end);
f1=f(1:n1);

%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %index of the nodes in the element
    %disp(Te);
    Xe=X(Te,:); %coordinates of the nodes in the element
    n=size(Xe,1); %nombre de neouds dans l'élément numéro e 
   
    [Ke,fe]=elementalComputations(Xe,IPcoord,IPweights,N,Nxi,Neta,sourceTermFunction);
  
    Ke3=Ke(n+1:end,1:n);
    Ke1=Ke(1:n,1:n);
    Ke4=Ke(n+1:end,n+1:end);
    Ke2=Ke(1:n,n+1:end);

    fe1=fe(1:n);
    fe2=fe(n+1:end);
  
    f1(Te)=f1(Te)+fe1;
    f2(Te)=f2(Te)+fe2;

    K1(Te,Te)=K1(Te,Te)+Ke1;
    K2(Te,Te)=K2(Te,Te)+Ke2;
    K3(Te,Te)=K3(Te,Te)+Ke3;
    K4(Te,Te)=K4(Te,Te)+Ke4;
end
K=[K1 K2;K3 K4]; f=[f1;f2];

%_______________________________________
%Calcul de la matriu i vector elementals
function [Ke,fe]=elementalComputations(Xe,IPcoord,IPweights,N,Nxi,Neta,sourceTermFunction)

E=2.5;
mu=0.25;

nnodes = size(Xe,1); %number of nodes in the element 
Ke=zeros(2*nnodes); %size of the elemental matrix
fe=zeros(2*nnodes,1);
xe = Xe(:,1); ye = Xe(:,2);

%Bucle en punts d integraci
for k=1:length(IPweights)
    Nk=N(k,:);  %k-th quadrature points on the element
    Nkxi=Nxi(k,:);
    Nketa=Neta(k,:); % derivatives with respect to eta in the reference element evaluated on the gauss point 
    xk = Nk*Xe; %map the gauss point to the physical element
    %Jacobian of the isoparametric change of variable 
    J = [Nkxi*xe Nkxi*ye;Nketa*xe Nketa*ye];
    % Derivadas de las funciones de forma respecto a (x,y)
    Nkxy = J\[Nkxi;Nketa]; %dérivative of the shape function with respect to x and y 
    Nkx=Nkxy(1,:);
    Nky=Nkxy(2,:);
    strain=[Nkx zeros(1,nnodes); zeros(1,nnodes) Nky; Nky Nkx]; %strain N 
    
    C=E/((1+mu)*(1-2*mu))*[1-mu mu 0 ; mu 1-mu 0; 0 0 1-2*mu];
    
    
    %diferencial de volum
    dxy=IPweights(k)*det(J);
    
    Ke = Ke + strain'*C*strain*dxy;%strain N^t*C*strain N 
    %fe = fe + sourceTermFunction(xk)*Nk'*dxy;
end
  