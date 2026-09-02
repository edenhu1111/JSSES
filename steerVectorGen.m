function a = steerVectorGen(sinTheta,sinPhi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:   steerVectorGen
% Description: 
% to generate steervector(UPA/ULA)
% Author: Eden HU (huyb@mail.sim.ac.cn)
% Input description
%
% sinTheta: the sine of azimuth angle
% sinPhi: the sine of vertical angle
% 
% Output description
% a: output steering function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Na Nv
nn1 = (0:Na-1).';
nn2 = (0:Nv-1).';
aa = exp(-1j*pi*sinTheta*nn1);
av = exp(-1j*pi*sinPhi  *nn2);
a = kron(aa,av);
end