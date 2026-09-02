function I = vec2imag(vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: vec2imag
% Description: transform the sparse vector into image
% Author: Eden Hu (huyb@mail.sim.ac.cn)
%
% Input description
% vec: the image vector
%
% Output description
% I: recovered image
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Nsize
I = reshape(vec,Nsize,Nsize);
end