function q=gradientSmoothing1D(u,X,T,theReferenceElement)

[nOfElements,nOfElementNodes]=size(T);
nOfNodes = size(X,1);
N=theReferenceElement.N; %basis functions at Gauss points
dNdxi=theReferenceElement.dNdxi;
wIP=theReferenceElement.IPweights';

M=spalloc(nOfNodes,nOfNodes,3*nOfNodes);
b=zeros(nOfNodes,1);
%Loop in elements
for e=1:nOfElements
    Te=T(e,:); %nodes of the element
    Xe=X(Te); %coordinates of the nodes in the element
    J=dNdxi*Xe; %Jacobian at Gauss points
    dNdx = diag(1./J)*dNdxi; %derivatives with respect to x
    dx=wIP.*J;
    M(Te,Te)=M(Te,Te) + N'*(diag(dx)*N); %assembly
    b(Te) = b(Te) + N'*(dx.*(dNdx*u(Te)));
end

q=M\b;
