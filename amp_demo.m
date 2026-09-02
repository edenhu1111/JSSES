clc;
clear all;
close all;
global gg HH
load Phi4test.mat
N = 400; % length of vector to be recovered
M = 160; % number of measurement

% phi = accuratePhi/(norm(accuratePhi,2)+0.01);
% [M,N] = size(phi);
% [U0,S,V0] = svd(phi);
% s = diag(S); s = s(s>=0.001);rank = length(s);
% S = diag(s);

% phirand = (1/sqrt(2))*(normrnd(0,1/sqrt(M),M,N) + 1i*normrnd(0,1/sqrt(M),M,N)); % Sensing matrix construction for theroetical bound
% [~,~,V] = svd(phirand);
% % V = dftmtx(N);
% % ind = randperm(N,rank);
% U0 = U0(:,1:rank); 
% V = V(:,1:rank);
% % ind = randperm(length(1:rank));
% % V = V(:,ind);
% 
% phi = U0*S*V';

phi = (1/sqrt(2))*(normrnd(0,1/sqrt(M),M,N) + 1i*normrnd(0,1/sqrt(M),M,N)); % Sensing matrix construction for theroetical bound
[U,~,V] = svd(phi);
rank = round(0.6*M);
S = diag(logspace(0,log10(1),rank));
phi = U(:,1:rank)*S*V(:,1:rank)';
phi = phi/norm(phi,2);
% 
% sigmaPhi = 2/sqrt(M*N);
% phiHat = phi + (1/sqrt(2)*normrnd(0,sigmaPhi,size(phi,1),size(phi,2)) + ...
%     1j/sqrt(2)*normrnd(0,sigmaPhi,size(phi,1),size(phi,2)));
HH = phi;
% 
% phi = dftmtx(N);
% ind = randperm(N,M);
% % ind = 1:5:N; M = length(ind);
% ind = sort(ind,'ascend');
% phi = phi(ind,:);
% phi = phi/norm(phi);
% [U,S,V] = svd(phi);
% rankP = round(0.56*M);
% phi = U(:,1:rankP)*S(1:rankP,1:rankP)*(V(:,1:rankP))';
% phi = phi/norm(phi);

% phi = normrnd(0,1/sqrt(M),M,N);
x = zeros(N,1); % Initializing sparse vector to be recovered
k = 10; % Sparsity level
uset = randperm(N,k); 
x(uset) = (rand(k,1) + 1i*rand(k,1))*1e0; % Sparse vector initialized
gg = x;
% x(uset) = (rand(k,1))*1e1; % Sparse vector initialized
% x = randn(N,1)*1e+1;
noise = sqrt(1/2)*(normrnd(0,1,M,1) + 1i*normrnd(0,1,M,1)); % zero mean, unit covariance complex noise vector
% noise = (normrnd(0,1,M,1));
var = 1e-5;
noise = sqrt(var)*noise;
y = phi*x + noise; % create measurement
niter = 50; % nu`mber of iteration
%% Approximate Message Passing for basis selection
% xest = amp(y,phi,niter);
% [xest] = vamp(y,phi,1/var,niter);
% [xest] = vampSVD(y,phi,1/var,niter);
close all;

inputParam.muPrior = 0;
inputParam.Lambda = 0.99;
inputParam.gammaOmega = 1/(norm(y,2)^2/(101*M));
% inputParam.gammaOmega = 1e+3;
inputParam.gammaPrior = norm(phi,2)^2*inputParam.Lambda/(norm(y,2)^2 - M/inputParam.gammaOmega);
inputParam.Normalization = 1;
inputParam.niter = 10;
inputParam.EMiter = 100;
inputParam.convBreaker = false;
inputParam.dampFac = 0.8;
inputParam.dampFacGam = 0.5;
[xest,state] = EMBGvampSVD(y,phi,inputParam);
NMSE = state.NMSE;
NMSEz = state.NMSEz;
% [xest] = laplavamp(y,phi,1/var,niter);

figure(1);
stem(abs(x),'x');hold on;stem(abs(xest),'o');
legend('Ground Truth','Est',Location='best');
% figure(2);semilogy(reshape(NMSE,[],1));xlabel('Iteration',Interpreter='latex');
figure(2);semilogy(NMSE(end,:));xlabel('Iteration',Interpreter='latex');

ylabel('NMSE',Interpreter='latex');
% figure(3);semilogy(reshape(NMSEz,[],1));xlabel('Iteration',Interpreter='latex');
% ylabel('NMSEz',Interpreter='latex');

% [abs(xest) abs(x)]

% %% Bi-Approximate Message Passing for basis selection
% % close all;
% xUAVEst = [];
% vecTOEst = [];
% xR = [];
% subcarrInd = [];
% sigmaPos = [];
% inputParam.sMDeriv.Phidx = [];
% inputParam.sMDeriv.Phidy = [];
% inputParam.sMDeriv.Phidz = [];
% inputParam.sMDeriv.PhiddelaywLoS = [];
% 
% 
% 
% inputParam.muPrior = 0;
% inputParam.Lambda = 0.99;
% inputParam.Normalization = 1;
% 
% inputParam.gammaOmega = 1/(norm(y,2)^2/(101*M));
% % inputParam.gammaOmega = 1e+3;
% inputParam.gammaPrior = norm(phi,2)^2*inputParam.Lambda/(norm(y,2)^2 - M/inputParam.gammaOmega);
% % inputParam.vH = 1/sigmaPhi;
% inputParam.vH = +inf;
% 
% inputParam.EMiter = 100;
% 
% inputParam.niter = 1;  % must be 1
% inputParam.niterA = 15;
% 
% inputParam.niterC = 1;
% inputParam.niterD = 10;
% 
% inputParam.convBreaker = false;
% inputParam.dampFac = 1;
% inputParam.dampFacGam = 1;
% inputParam.xUAV = xUAVEst;
% inputParam.Tdelay = vecTOEst;
% inputParam.xR = xR;
% inputParam.subcarrInd = subcarrInd;
% inputParam.stepSize = 1;
% inputParam.sigmaPos = sigmaPos;
% 
% inputParam.armijoBeta = 0.2;
% inputParam.armijoSigma = 0.01;
% 
% [xest2,state2] = BiTMP3(y,phiHat,inputParam);
% NMSE = state2.NMSE;
% % [xest] = laplavamp(y,phi,1/var,niter);
% 
% figure(3);
% stem(abs(x),'x');hold on;stem(abs(xest2),'o');hold off;
% legend('Ground Truth','Est',Location='best');
% % figure(2);semilogy(reshape(NMSE,[],1));xlabel('Iteration',Interpreter='latex');
% figure(4);semilogy(NMSE(end,:));xlabel('Iteration',Interpreter='latex');
% ylabel('NMSE',Interpreter='latex');title('NMSE',Interpreter='latex');
% 
% figure(5);semilogy(reshape(state2.NMSEH,[],1));xlabel('Iteration',Interpreter='latex');
% ylabel('NMSEH',Interpreter='latex');title('NMSEH',Interpreter='latex');ylim([1e-2 2e0]);
% 
% figure(6);semilogy(state2.gammaOmega);xlabel('Iteration',Interpreter='latex');
% ylabel('$\gamma_\omega$',Interpreter='latex');title('Variance',Interpreter='latex');
% % [abs(xest) abs(x)]
