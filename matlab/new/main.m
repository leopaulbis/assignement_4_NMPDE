%% SOLUTION OF THE LAPLACE EQUATION IN 1D: -u_xx=f with Dirichlet boundary conditions
close all; clear all;

degree = 1;
%Preprocess (mesh)
n=30; %number of elements
[X,T]=createUniformMesh1D(0,1,n,degree); %nodal coordinates and connetivity matrix

%Reference element
theReferenceElement = referenceElement1D(degree);

%System computation
[K,f]=computeSystem(X,T,@sourceTerm,theReferenceElement);


%Dirichlet boundary conditions (system reduction)
u0=0; u1=0;
f1=K(2:end-1,1); fend=K(2:end-1,end);
f=f(2:end-1)-u0*f1-u1*fend;
K=K(2:end-1,2:end-1);

%System solution
u=K\f;

%Postprocess
u=[u0;u;u1];
im=find(abs(X-0.5)<1.e-5); if ~isempty(im),fprintf('   u^h(0)=%f\n',u(im)); end
figure(1), plotSolution1D(X,T,u), title('u^h')
figure(2), clf, plotDerivative1D(X,T,u), title('derivative')
 
%Error estimate
%Gradient smoothing: L2 projection of the gradient
q=gradientSmoothing1D(u,X,T,theReferenceElement);
figure(3), hold on, plotSolution1D(X,T,q), hold off



