function xest = amp(y,phi,niter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: amp
% Description: AMP for solving linear function 
% Author: Eden Hu (huyb@mail.sim.ac.cn)
%
% Input description
% y: the noisy data
% phi: sensing matrix
% niter: maximum iteration number
%
% Output description
% xest: estimated x
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(phi,2);
M = size(phi,1);
%% Approximate Message Passing for basis selection
% Initializing
xest = zeros(N,1);
r = randn(N,1) + 1j*randn(N,1);
gamma = 100;
v = zeros(M,1);
eta = @(x,beta) (x./abs(x)).*(abs(x)-beta).*(abs(x)-beta > 0); % denoising function
dEta = @(x,beta) (abs(x)-beta > 0);
a = 0.7;
for iter=1:niter
   xest = eta(r, gamma);
   alpha = mean(dEta(r,gamma));
   v = y - phi*xest + N/M*alpha*v;
   r = xest + phi'*v;
   gamma = a*gamma +(1-a) * M/(norm(v,'fro'))^2;
end

%[abs(xest) abs(x)]
end