function [K,f]=computeSystem(X,T,sourceTerm,theReferenceElement)

[nOfElements,nOfElementNodes]=size(T);
nOfNodes = size(X,1);

K=spalloc(nOfNodes,nOfNodes,3*nOfNodes);
f=zeros(nOfNodes,1);
%Loop in elements
for e=1:nOfElements
    Te=T(e,:); %nodes of the element
    Xe=X(Te); %coordinates of the nodes in the element
    [Ke,fe]=elementalSourceTerm(Xe,sourceTerm,theReferenceElement);
    K(Te,Te)=K(Te,Te) + Ke; %assembly
    f(Te) = f(Te) + fe;
end


%__Computation of the source term with numerical integration
function [Ke,fe]=elementalSourceTerm(Xe,sourceTerm,theReferenceElement)

N=theReferenceElement.N; %basis functions at Gauss points
dNdxi=theReferenceElement.dNdxi;
J=dNdxi*Xe; %Jacobian at Gauss points
dNdx = diag(1./J)*dNdxi; %derivatives with respect to x
wIP=theReferenceElement.IPweights';
dx=wIP.*J;
xIP= N*Xe; %x-coordinates of the integration points
fIP=sourceTerm(xIP); %source term at Gauss points
fe=N'*(dx.*fIP);
Ke=dNdx'*(diag(dx)*dNdx); 
