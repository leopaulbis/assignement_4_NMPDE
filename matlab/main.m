%PREPROCES
%Reference element
clear all;
degree = 1; typeOfElement=1; %1=TRI, 0=QUA
theReferenceElement = createReferenceElement(degree,typeOfElement);
nOfElementNodes = size(theReferenceElement.N,2);
%figure(1), drawReferenceElement(theReferenceElement);
%Mesh: regular mesh in a rectangular domain [0,1]x[0,1]
[X,T] = CreateMesh(1,nOfElementNodes,[0,1,0,1],15,15);
%disp(X);
figure(2), clf
PlotMesh(T,X,typeOfElement,'k-',1);
%Definition of Dirichlet boundary conditions
x = X(:,1); y = X(:,2); tol=1.e-10;
nodesCCD = find(abs(x)<tol|abs(x-1)<tol|abs(y)<tol|abs(y-1)<tol);%index of the nodes on the boundary
%hold on, plot(x(nodesCCD),y(nodesCCD),'bo','MarkerSize',16); hold off

uCCD=DirichletValue(X(nodesCCD,:));
uCCD=[uCCD;zeros(size(uCCD,1),1)]; %value 0 for the y component and other the x components. 

%System of equations (without BC)
[K,f]=computeSystemLaplace(X,T,theReferenceElement,@sourceTerm);
n=size(X,1);

%Boundary conditions & SYSTEM SOLUTION
methodDBC = 1; % 1=System reduction or 2=Lagrange multipliers 
if methodDBC == 1 %System reduction
    nodesCCD=[nodesCCD ;n+nodesCCD];
    unknowns= setdiff([1:n+size(X,1)],nodesCCD);%indice des noeuds qui ne sont pas sur le bord de dirichlet 
    
    f = f(unknowns)-K(unknowns,nodesCCD)*uCCD;
    K=K(unknowns,unknowns);
    %System solution
    sol=K\f;
    %Nodal values: system solution and Dirichlet values
    u = zeros(2*size(X,1),1);
    u(unknowns) = sol; u(nodesCCD) = uCCD;
else %LagrangeMultipliers
    nOfDirichletDOF = length(uCCD); nOfDOF = size(K,1);
     nodesCCD=[nodesCCD ;n+nodesCCD];
    A = spalloc(nOfDirichletDOF,nOfDOF,nOfDirichletDOF);
    A(:,nodesCCD) = eye(nOfDirichletDOF);
    b = uCCD;
    Ktot = [K A'; A spalloc(nOfDirichletDOF,nOfDirichletDOF,0)];
    ftot = [f;b];
    disp(ftot);
    sol = Ktot\ftot;
    u = sol(1:nOfDOF); lambda = sol(nOfDOF+1:end);
end
%disp(u);
u_1=u(1:n);
u_2=u(n+1:end);
%POSTPROCESS
% figure(3)
% PlotNodalField(u_2,X,T), title('FEM solution')

% figure(4)
% PlotNodalField(u_2,X,T), title('FEM solution')

% X=X+[u_1,u_2];
% figure;
% PlotMesh(T,X,typeOfElement,'k-',1);

E=2.5;
mu=0.25;
C=E/((1+mu)*(1-2*mu))*[1-mu mu 0 ; mu 1-mu 0; 0 0 1-2*mu];

[du_1dx,du_1dy]=computeGradientSmoothing(u_1,X,T,theReferenceElement);
[du_2dx,du_2dy]=computeGradientSmoothing(u_2,X,T,theReferenceElement);

epsilon=[du_1dx';du_2dy';du_1dy'+du_2dx'];
sigma=zeros(3,size(u_1,1));

for i=1:size(epsilon,2)
    sigma(:,i)=C*epsilon(:,i);
end 

sigma_VM=zeros(1,size(u_1,1));

for i=1:size(u_1,1)
    sigma_VM(1,i)=sqrt(3/2*(sigma(1,i)^2+sigma(2,i)^2-sigma(1,i)*sigma(2,i)+3*sigma(3,i)^2));
end 

X=X+[u_1,u_2];
figure;
PlotMesh(T,X,typeOfElement,'k-',1);
%disp(sigma_VM);
figure;
PlotNodalField(sigma_VM',X,T);

%view(2)

% sigma=epsilon;
% figure(4);
% 
% plotSolution1D(X,T,sigma);






