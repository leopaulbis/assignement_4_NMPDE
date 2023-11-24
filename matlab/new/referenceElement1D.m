function theReferenceElement = referenceElement1D(degree)
%
% theReferenceElement = referenceElement1D(degree)
%

switch degree
    case 1
        GaussCoordinates = [-sqrt(3)/3; sqrt(3)/3];
        GaussWeights = [1 1];
        N = [(1-GaussCoordinates)/2 (1+GaussCoordinates)/2];
        dNdxi = [-1/2 1/2; -1/2 1/2];
    case 2
        GaussCoordinates = [-sqrt(3/5); 0; sqrt(3/5)];
        GaussWeights = [5 8 5]/9;
        N = [(GaussCoordinates-1).*GaussCoordinates/2 (1+GaussCoordinates).*GaussCoordinates/2 1-GaussCoordinates.^2];
        dNdxi = [GaussCoordinates-1/2 GaussCoordinates+1/2 -2*GaussCoordinates];    
    otherwise
        error(['No se puede utilizar una interpolacion de grado ',num2str(degree)]);
end

theReferenceElement.IPcoord=GaussCoordinates;
theReferenceElement.IPweights=GaussWeights;
theReferenceElement.N = N;
theReferenceElement.dNdxi = dNdxi;

