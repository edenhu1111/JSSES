clc;
clear all;
close all;
global gg HH
% N = 421; % length of vector to be recovered
% M = 4000; % number of measurement
load('dataset\3\Phid.mat');
load('dataset\3\Phip.mat');
load('dataset\3\PhidGenie.mat');
load('dataset\3\PhipGenie.mat');
load("dataset\3\d.mat");
load("dataset\3\noisedatap.mat");
load("dataset\3\noisedatad.mat");
load("dataset\xOG.mat");
load("dataset\3\zp.mat");
load("dataset\3\zd.mat");
% phi = repmat(kron(ones(32,1),d(:,1,3)),1,size(Phid,2)).*Phid(:,:,1);
% phii = Phip(:,:,1);

phi = repmat(kron(ones(16,1),d(:,1,3)),1,size(Phid,2)).*PhidGenie(:,:,1);
phii = Phip(:,:,1);

y  = zd(:,1,3) + 0.97e-2*linspace(1,1,size(Phid,1))'.*noisedatad(:,1,1);
yy = zp(:,1)   + 0.97e-2*linspace(1,1,size(Phip,1))'.*noisedatap(:,1,1);
% phi = 1/sqrt(2)*(randn(M,N) + 1j*randn(M,N)); % Sensing matrix construction for theroetical bound
% % [U,~,V] = svd(phi);rank = round(N);
% % S = diag(logspace(0,log10(0.01),rank));
% % phi = U(:,1:rank)*S*V(:,1:rank)';
% phi = phi/norm(phi,'fro');
% scale = 10;


% x = zeros(N,1); % Initializing sparse vector to be recovered
% k = 10; % Sparsity level
% uset = randperm(N,k); 
% x(uset) = (rand(k,1) + 1i*rand(k,1))*1e0; % Sparse vector initialized
gg = x(:,1);
% x(uset) = (rand(k,1))*1e1; % Sparse vector initialized
% x = randn(N,1)*1e+1;
% noise  = sqrt(1/2)*(normrnd(0,1,M,1)       + 1i*normrnd(0,1,M,1)); % zero mean, unit covariance complex noise vector
% noisee = sqrt(1/2)*(normrnd(0,1,M/scale,1) + 1i*normrnd(0,1,M/scale,1)); % zero mean, unit covariance complex noise vector
% % noise = (normrnd(0,1,M,1));
% var = 1e-6;
% noise = sqrt(var)*noise;
% y = phi*x + noise; % create measurement
% phi0 = phi;
% phi(:,uset) = phi(:,uset) + 0.00003*(repmat(linspace(3,16,M)',1,k)).*(randn(size(phi0,1),k) + 1j*randn(size(phi0,1),k));
% yy = y(1:scale:M);
% phii = phi(1:scale:M,:);

niter = 50; % number of iteration
%% SBL for basis selection

inputParam.gammaA = 1e-8;
inputParam.gammaB = 1e-8;
inputParam.gammaC = 1e-8;
inputParam.gammaD = 1e-8;
inputParam.rho0 = 1e4;
inputParam.Niter = 100;
eyeS = eye(32);

inputParam.constell = 0;
inputParam.eyeS = eyeS;
inputParam.convBreak = false;
inputParam.gammaOmega = 1e3;
inputParam.Niter = 10;
inputParam.NiterA = 20;
inputParam.lambda0 = 0.5;
inputParam.dataNum = 10;
inputParam.JESD = 0;
inputParam.lambdaS = 0.001;
[xest1,~,state1] = VBI2(yy,y,phii,phi,inputParam);
NMSE1 = vec(state1.NMSE);
NMSE1 = NMSE1(NMSE1~=0);
% [xest1,state1] = sbl(y,phi,inputParam);
% NMSE1 = state1.NMSE;

inputParam.gammaA = 1e-8;inputParam.gammaB = 1e-8;

[xest2,state2] = sbl(yy,phii,inputParam);
NMSE2 = state2.NMSE;


inputParam.SV = 0;  
inputParam.convBreak = false;
inputParam.gammaOmega = 1e3;
inputParam.Niter = 10;
inputParam.NiterA = 20;
inputParam.lambda0 = 0.001;
[xest2,~,state2] = SVVBI(yy,y,phii,phi,inputParam);
NMSE2 = vec(state2.NMSE);
NMSE2 = NMSE2(NMSE2~=0);
[xest3,~,state3] = SVVBI(y,yy,phi,phii,inputParam);

NMSE3 = vec(state3.NMSE);
NMSE3 = NMSE3(NMSE3~=0);
%%
figure(1);
subplot(1,3,1);
stem(abs(gg),'x');hold on;stem(abs(xest1),'o');
legend('Ground Truth','Est',Location='best');
% 
subplot(1,3,2);
stem(abs(gg),'x');hold on;stem(abs(xest2),'o');
legend('Ground Truth','Est',Location='best');

subplot(1,3,3);
stem(abs(gg),'x');hold on;stem(abs(xest3),'o');
legend('Ground Truth','Est',Location='best');
%%
figure(3);
semilogy(NMSE1);
hold on;
semilogy(NMSE2);
semilogy(NMSE3);legend('1','2','3');
xlabel('Iteration',Interpreter='tex');
ylabel('NMSE',Interpreter='tex');
% figure(4);semilogy(abs(state1.rho));hold on;semilogy(abs(state2.rho));

