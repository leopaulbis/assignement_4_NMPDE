function [ x,w ] = getTriangleQuadrature( N )

    [X,Y,Wx,Wy]=triquad(N,[0 0 ; 1 0 ; 0 1]);
    
    x = zeros(N*N,2);
    w = zeros(N*N,1);
    
    count = 0;
    
    for i=1:N
        for j=1:N
            count=count+1;
            x(count,:)=[X(i,j) Y(i,j)];
            w(count)=Wx(i)*Wy(j);
        end
    end  
    
end

