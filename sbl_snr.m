clc;
clear all;
close all;
global gg HH
load Phi4test.mat
N = 200; % length of vector to be recovered
M = 400; % number of measurement



phi = (1/sqrt(2))*(normrnd(0,1/sqrt(M),M,N) + 1i*normrnd(0,1/sqrt(M),M,N)); % Sensing matrix construction for theroetical bound
[U,~,V] = svd(phi);rank = round(0.3*N);
S = diag(logspace(0,log10(0.0001),rank));
phi = U(:,1:rank)*S*V(:,1:rank)';
phi = phi;

% 

x = zeros(N,1); % Initializing sparse vector to be recovered
k = 10; % Sparsity level
uset = randperm(N,k); 
x(uset) = (rand(k,1) + 1i*rand(k,1))*1e0; % Sparse vector initialized
gg = x;
% x(uset) = (rand(k,1))*1e1; % Sparse vector initialized
% x = randn(N,1)*1e+1;
noise = sqrt(1/2)*(normrnd(0,1,M,1) + 1i*normrnd(0,1,M,1)); % zero mean, unit covariance complex noise vector
% noise = (normrnd(0,1,M,1));
VAR = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5];
VARdB = 10*log10(VAR);
% var = 1e-4;
% noise = sqrt(var)*noise;
for ii = 1:length(VAR)
    var = VAR(ii);
    noise = sqrt(var)*noise;

%% SBL for basis selection
    y = phi*x + noise; % create measurement
    niter = 50; % number of iteration
    inputParam.gammaA = 1e-8;
    inputParam.gammaB = 1e-8;
    inputParam.gammaC = 1e0;
    inputParam.gammaD = 1e-8;
    inputParam.rho0 = 1e4;
    inputParam.Niter = 100;
    
    
    
    [xest1,~] = sbl(y,phi,inputParam);
    NMSE1(ii) = norm(xest1 - x,2)^2/norm(x,2)^2;
    
    inputParam.gammaA = 1;inputParam.gammaB = 1e-8;
    
    [xest2,state2] = sbl(y,phi,inputParam);
    NMSE2(ii) = norm(xest2 - x,2)^2/norm(x,2)^2;
    
    inputParam.convBreak = false;
    inputParam.gammaOmega = 1e4;
    inputParam.Niter = 5;
    inputParam.NiterA = 20;
    inputParam.lambda0 = 0.02;
    [xest3,state3] = VBI(y,phi,inputParam);
    NMSE3(ii) = norm(xest3 - x,2)^2/norm(x,2)^2;
%     NMSE3 = NMSE3(NMSE3~=0);
end
%%
figure(1);
semilogy(VARdB,NMSE1,'-o');hold on;
semilogy(VARdB,NMSE2,'--s');semilogy(VARdB,NMSE3,':<');
legend('1','2','3');grid on;

