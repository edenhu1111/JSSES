function [output] = measurmentMatrixGen(input)
grid = input.grid;
M = input.M; N = input.N;
Phi = zeros(M,N);
ind = (0:M-1).';
for nn = 1:N
    Phi(:,nn) = exp(1j*grid(nn)*ind);
    Phi1(:,nn) = 1j*ind.*   Phi(:,nn);
    Phi2(:,nn) = -1*ind.^2.*Phi(:,nn);
end
sf = sqrt(M);
output.Phi = Phi/sf;
output.Phi1 = Phi1/sf;
output.Phi2 = Phi2/sf;

end