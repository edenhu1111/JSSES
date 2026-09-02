function [Phi,Phi1,Phi2] = sMGen(dim1,dim2,grid)
    gridLen = length(grid(1,:));
    for ii = 1:gridLen
        aa = exp(-1j*pi*(0:dim1-1)'*grid(1,ii)); bb = exp(-1j*pi*(0:dim2-1)'*grid(2,ii));
        Phi(:,ii) = kron(aa,bb);
        a1 = -1j*pi*(0:dim1-1)'.*aa; b1 = -1j*pi*(0:dim2-1)'.*bb;
        
        Phi1(:,ii) = kron(a1,bb);
        Phi2(:,ii) = kron(aa,b1);
    end
    scalingFac = norm(Phi,'fro')/sqrt(size(grid,2));
    Phi = Phi/scalingFac;
    Phi1 = Phi1/scalingFac;
    Phi2 = Phi2/scalingFac;
end