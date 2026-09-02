%% Initalize the environment
clc;
clear;close all;
dbstop if error;

global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle gg xRange HH
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 30e+9;                    % carrier frequency (Ka)

deltaF = 120e+3;   % subcarrier spacing
BW = deltaF*1024;  
T = 1 / deltaF;   % symbol duration
Tcp = T / 4;      % cyclic prefix duration
Ts = T + Tcp;     % total symbol duration
sigmaPos = 1/sqrt(3);
sigma2 = 10^(-65/10);
xRange = 50;
niter = 20;
Na = 64;
Nv = 1;
Nsize = 20;
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
xR = [-50; -0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = size(xR,2);
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
subcarrInd = 0:64:1023;
% subcarrInd = sort(randperm(1024,8));
% pilot = zadoffChuSeq(17,length(subcarrInd));
pilot = ones(length(subcarrInd),1);
xUAV = [50;0;60];
Nt = size(xUAV,2);
% xUAV = [-10,  0, 10,-10, 0 , 10,-10, 0, 10;...
%         -10,-10,-10,  0, 0 , 0 , 10, 10,10;...
%          60, 60, 60, 60, 60, 60, 60, 60,60];
% xUAV = [ -5,  0,  5;...
%           0,  0,  0;...
%          60, 60, 60]; 
% xUAV = [-40,  0, 40,-40, 0 , 40,-40, 0, 40;...
%         -40,-40,-40,  0, 0 , 0 , 40, 40,40;...
%          60, 60, 60, 60, 60, 60, 60, 60,60];
% xUAV = xUAV + 1*randn(size(xUAV));
xUAVEst = xUAV + sigmaPos*(rand(size(xUAV))-1/2);
% xUAVEst = xUAV;

% vecTO = 0*1/BW*(rand(Nr,1) - 1/2)/10;
vecTO = 1/BW*(rand-1/2)*ones(Nr,1)/2;
% vecTO = 1/BW*(rand-1/2);

dArray = c0/fc/2;
biasAngle = [0,0,0,0]';                                                  %degree
% biasAngle = [90]';                                                  %degree

% biasAngle = zeros(Nr,2); 
% vecTO = zeros(Nr,1);
vecTOEst = zeros(size(vecTO));
% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
xTarget = zeros(3,16);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation
% g = 100* (randn(size(xTarget,2)) + 1j*randn(size(xTarget,2)));
% accuratePhi = sensingMatrixGen(xUAV,xR,subcarrInd,vecTO);
accuratePhi = sensingMatrixGenWPilot(xUAV,xR,subcarrInd,pilot,vecTO);
% accuratePhi = accuratePhi./(ones(size(accuratePhi,1),1)*vecnorm(accuratePhi,1));
% accuratePhi = accuratePhi/norm(accuratePhi,2);
% accuratePhi = normrnd(0,1/sqrt(2),size(accuratePhi))+1j*normrnd(0,1/sqrt(2),size(accuratePhi));
% accuratePhi = accuratePhi/norm(accuratePhi,2);
HH = accuratePhi;
% inaccuratePhi = sensingMatrixGen(xUAVEst,xR,subcarrInd,vecTOEst);
[inaccuratePhi,sensingMatrixDeriv]...
    = sensingMatrixGenWPilot(xUAVEst,xR,subcarrInd,pilot,vecTOEst);
% inaccuratePhi =
% inaccuratePhi./(ones(size(inaccuratePhi,1),1)*vecnorm(inaccuratePhi,2));


% nrm = norm(inaccuratePhi,2);
% inaccuratePhi  = inaccuratePhi/nrm;
% 
% sensingMatrixDeriv.Phidx = sensingMatrixDeriv.Phidx/nrm;
% sensingMatrixDeriv.Phidy = sensingMatrixDeriv.Phidy/nrm;
% sensingMatrixDeriv.Phidz = sensingMatrixDeriv.Phidz/nrm;
% sensingMatrixDeriv.PhiddelaywLoS = sensingMatrixDeriv.PhiddelaywLoS/nrm;


% % sensingMatrixDeriv.Phidxdx  = sensingMatrixDeriv.Phidxdx/norm(inaccuratePhi,2);
% % sensingMatrixDeriv.Phidxdx2 = sensingMatrixDeriv.Phidxdx2/norm(inaccuratePhi,2);
% % sensingMatrixDeriv.Phidddd  = sensingMatrixDeriv.Phidddd/norm(inaccuratePhi,2);
% % sensingMatrixDeriv.Phidddd2 = sensingMatrixDeriv.Phidddd2/norm(inaccuratePhi,2);



% 
% g = zeros(size(accuratePhi,2),1);
% g(19900) = 100;
% r = accuratePhi*g + sqrt(sigma2)*(randn(size(accuratePhi,1),1) + 1j*randn(size(accuratePhi,1),1));
% g = 3*(randn(size(xTarget,2),1) + 1j*randn(size(xTarget,2),1));
g = normrnd(1/sqrt(2),0.1,size(xTarget,2),1) + 1j*normrnd(1/sqrt(2),0.1,size(xTarget,2),1);
g = exp(1j*2*pi*randn(size(g))).*g;

originalImage = zeros(Nsize,Nsize);
for ii = 1:4
    for jj = 1:4
        originalImage(round(Nsize/5*ii),round(Nsize/5*jj)) = g((ii-1)*4+jj);
    end
end

gg = reshape(originalImage,[],1);

% gg = zeros(Nsize*Nsize,1);
% ind = randperm(Nsize*Nsize,size(xTarget,2));
% gg(ind) = g;

gg = [10*ones(Nr,1);gg];
% gg = [10*ones(Nr,1);zeros(size(gg))];
r = accuratePhi*gg + (normrnd(0,sqrt(sigma2/2),size(accuratePhi,1),1) +...
    1j*normrnd(0,sqrt(sigma2/2),size(accuratePhi,1),1));
% r = accuratePhi*gg;
% [r] = rxSignalGen(xUAV,xR,xTarget,g,subcarrInd,vecTO);
% [r] = rxSignalGenWPilot(xUAV,xR,xTarget,g,subcarrInd,pilot,vecTO);

% SNR = 10*log10(norm(accuratePhi*gg,2)^2/norm(r-accuratePhi*gg,2)^2)

% accuratePhiEq = [real(accuratePhi),-imag(accuratePhi);imag(accuratePhi),real(accuratePhi)];
% inaccuratePhiEq = [real(inaccuratePhi),-imag(inaccuratePhi);imag(inaccuratePhi),real(inaccuratePhi)];
% rEq = [real(r);imag(r)];
% gg = [real(gg);imag(gg)];
[M,N] = size(accuratePhi);
%% BG-VAMP
% xest1 = amp(rEq,accuratePhiEq,1400,niter);
% xest1 = xest1(1:length(xest1)/2)+1j*xest1(length(xest1)/2+1:end);
% xest1 = amp(r,accuratePhi,niter);
inputParam.muPrior = [0*ones(Nr,1);zeros(size(accuratePhi,2)-Nr,1)];
inputParam.Lambda = 0.999;

% inputParam.muPrior = 0;
inputParam.gammaOmega = 1/(norm(r,2)^2/(101*M));
inputParam.gammaPrior = norm(accuratePhi,2)^2*inputParam.Lambda/(norm(r,2)^2 - M/inputParam.gammaOmega);
inputParam.Normalization = 1;
inputParam.niter = 10;
inputParam.EMiter = 50;
inputParam.convBreaker = false;
inputParam.dampFac = 0.6;
inputParam.dampFacGam = 0.5;
% xest1 = vampSVD(r,  accuratePhi,niter,inputParam);
% [xest1,state] = EMBGvampSVD(r,accuratePhi,inputParam);
[xest1Origin,state1] = EMBGvampSVD(r,accuratePhi,inputParam);
xest1 = xest1Origin(Nt*Nr+1:end,:);
% [xest1,state] = EMBGvamp(r,accuratePhi,inputParam);
I1 = vec2imag(xest1);
figure(1);
subplot(1,2,1);
surf(abs(I1)/max(abs(I1),[],'all'));shading interp;
% [xest1,state] = EMBGvampSVD(rEq,accuratePhiEq,niter,gammaPrior,muPrior,Lambda,1/sigma2);
% xest1 = xest1(1:length(xest1)/2)+1j*xest1(length(xest1)/2+1:end);
subplot(1,2,2);
nmse1 = state1.NMSE;
spa = state1.sparsity;
nmsez = state1.NMSEz;
% semilogy(reshape(nmse1,[],1));
switch inputParam.convBreaker
    case false
        semilogy(reshape(nmse1(end,:),[],1));
    case true
        semilogy(reshape(nmse1,[],1));
end

%%
% xest1 = EMBGvamp(r,accuratePhi,niter,gammaPrior,muPrior,Lambda,1/sigma2);
inputParam.muPrior = [10*ones(Nr,1);zeros(size(accuratePhi,2)-Nr,1)];
inputParam.Lambda = 0.999;

% inputParam.muPrior = 0;
inputParam.gammaOmega = 1/(norm(r,2)^2/(101*M));
% inputParam.gammaOmega = 1e+2;
inputParam.gammaPrior = norm(accuratePhi,2)^2*inputParam.Lambda/(norm(r,2)^2 - M/inputParam.gammaOmega);
inputParam.Normalization = 1;

inputParam.niter = 10;
inputParam.EMiter = 50;
inputParam.convBreaker = false;
inputParam.dampFac = 0.5;
inputParam.dampFacGam = 0.5;
[xest2Origin,state2] = EMBGvampSVD(r,inaccuratePhi,inputParam);
xest2 = xest2Origin(Nt*Nr+1:end,:);


I2 = vec2imag(xest2);
% % % xest2 = amp(rEq,inaccuratePhiEq,1500,niter);
% % % xest2 = xest2(1:length(xest2)/2)+1j*xest2(length(xest2)/2+1:end);
% % % xest2 = amp(r,inaccuratePhi,niter);
% % % xest2 = vampSVD(r,inaccuratePhi,niter,Lambda,muPrior,gammaPrior,1/sigma2/1000^2);
% [xest2,nmse2] = EMBGvampSVD(r,inaccuratePhi,niter,gammaPrior,muPrior,Lambda,1/sigma2);
% % xest2 = EMBGvamp(r,inaccuratePhi,niter,gammaPrior,muPrior,Lambda,1/sigma2);
% 
% I2 = vec2imag(xest2);


% subplot(1,2,2);
figure(2);
subplot(1,2,1);
surf(abs(I2)/max(abs(I2),[],'all'));shading interp; 
subplot(1,2,2);
switch inputParam.convBreaker
    case false
        semilogy(reshape(state2.NMSEz(end,:),[],1));
    case true
        semilogy(reshape(state2.NMSEz,[],1));
end

%% EM-BiVAMP
 inputParam.muPrior = [10*ones(Nr,1);zeros(size(accuratePhi,2)-Nr,1)];
inputParam.Lambda = 0.999;
inputParam.gammaOmega = 1/(norm(r,2)^2/(101*M));
inputParam.gammaPrior = norm(accuratePhi,2)^2*inputParam.Lambda/(norm(r,2)^2 - M/inputParam.gammaOmega);
inputParam.Normalization = 1;

inputParam.EMiter = 50;
inputParam.EMOutiter = 50;
inputParam.niter = 10;  
inputParam.armijoBeta = 0.2;
inputParam.armijoSigma = 0.2;
  
inputParam.convBreaker = false;
inputParam.dampFac = 0.5;
inputParam.dampFacGam = 0.5;

inputParam.xUAV = xUAVEst;
inputParam.Tdelay = vecTOEst;

% inputParam.xUAV = xUAV;
% inputParam.Tdelay = vecTO;

inputParam.xR = xR;
inputParam.subcarrInd = subcarrInd;
inputParam.stepSize = 1e-12;
inputParam.sigmaPos = sigmaPos;

inputParam.sMDeriv = sensingMatrixDeriv;
tic
[xestEMOriginal,stateEM] = EMBivamp3(r,inaccuratePhi,inputParam);
% [xestEM,stateEM] = EMBivamp3(r,accuratePhi,inputParam);

% [xestEM,stateEM] = EMBGvampSVD(r,inaccuratePhi,inputParam);
toc
xestEM = xestEMOriginal(Nt*Nr+1:end,:);

% [xest1,state] = EMBGvamp(r,accuratePhi,inputParam);
I4 = vec2imag(xestEM);
figure(5);
subplot(1,2,1);
surf(abs(I4)/max(abs(I4),[],'all'));shading interp;
subplot(1,2,2);
switch inputParam.convBreaker
    case false
        semilogy(reshape(stateEM.NMSEz(end,end,:),[],1));
    case true
        semilogy(reshape(stateEM.NMSEz,[],1));
end

title('Bi-VAMP');

%%
SNR = 10*log10(norm(accuratePhi*gg,2)^2/norm(r-accuratePhi*gg,2)^2)
%%
%% EM-BiVAMP2
% gg = gg(1:Nt*Nr);
%  inputParam.muPrior = [0*ones(Nr,1)];
% inputParam.Lambda = 0.999;
% inputParam.gammaOmega = 1/(norm(r,2)^2/(11*M));
% inputParam.gammaPrior = norm(accuratePhi,2)^2*inputParam.Lambda/(norm(r,2)^2 - M/inputParam.gammaOmega);
% inputParam.Normalization = 0;
% 
% inputParam.EMiter = 50;
% inputParam.EMOutiter = 30;
% inputParam.niter = 10;  
% inputParam.armijoBeta = 0.2;
% inputParam.armijoSigma = 0.05;
% 
% inputParam.convBreaker = false;
% inputParam.dampFac = 0.6;
% inputParam.dampFacGam = 0.5;
% 
% inputParam.xUAV = xUAVEst;
% inputParam.Tdelay = vecTOEst;
% 
% % inputParam.xUAV = xUAV;
% % inputParam.Tdelay = vecTO;
% 
% inputParam.xR = xR;
% inputParam.subcarrInd = subcarrInd;
% inputParam.stepSize = 1e-12;
% inputParam.sigmaPos = sigmaPos;
% 
% inputParam.sMDeriv = sensingMatrixDeriv;
% tic
% [xestEM,stateEM] = EMBivamptest(r,inaccuratePhi(:,1:Nt*Nr),inputParam);
% % [xestEM,stateEM] = EMBivamp3(r,accuratePhi,inputParam);
% 
% % [xestEM,stateEM] = EMBGvampSVD(r,inaccuratePhi,inputParam);
% toc
% 
% 
% % [xest1,state] = EMBGvamp(r,accuratePhi,inputParam);
% figure(5);
% subplot(1,2,1);
% stem(abs(xestEM));hold on;stem(abs(gg));
% subplot(1,2,2);
% switch inputParam.convBreaker
%     case false
%         semilogy(reshape(stateEM.NMSEz(end,end,:),[],1));
%     case true
%         semilogy(reshape(stateEM.NMSEz,[],1));
% end
% 
% title('Bi-VAMP');


%%
% %% Proposed Algorithm
% inputParam.muPrior = [10*ones(Nr,1);zeros(size(accuratePhi,2)-Nr,1)];
% inputParam.Lambda = 0.999;
% inputParam.gammaOmega = 1/(norm(r,2)^2/(101*M));
% inputParam.gammaPrior = norm(accuratePhi,2)^2*inputParam.Lambda/(norm(r,2)^2 - M/inputParam.gammaOmega);
% inputParam.Normalization = 0;
% 
% inputParam.EMiter = 50;
% 
% inputParam.niter = 1;  % must be 1
% inputParam.niterA = 10;
% 
% inputParam.niterC = 1;
% inputParam.niterD = 1;
% 
% inputParam.convBreaker = false;
% inputParam.dampFac = 0.8;
% inputParam.dampFacGam = 0.5;
% inputParam.xUAV = xUAVEst;
% inputParam.Tdelay = vecTOEst;
% inputParam.xR = xR;
% inputParam.subcarrInd = subcarrInd;
% inputParam.stepSize = 1;
% inputParam.sigmaPos = sigmaPos;
% 
% inputParam.armijoBeta = 0.4;
% inputParam.armijoSigma = 0.4;
% 
% inputParam.sMDeriv = sensingMatrixDeriv;
% % xest1 = vampSVD(r,  accuratePhi,niter,inputParam);
% % [xest1,state] = EMBGvampSVD(r,accuratePhi,inputParam);
% tic
% [xestProposed,stateP] = BiTMP2(r,inaccuratePhi,inputParam);
% toc
% % [xestProposed,state1] = BiTMP(r,accuratePhi,inputParam);
% 
% xestProposed = xestProposed(Nt*Nr+1:end,:);
% % [xest1,state] = EMBGvamp(r,accuratePhi,inputParam);
% I3 = vec2imag(xestProposed);
% figure(3);
% % subplot(1,2,1);
% surf(abs(I3)/max(abs(I3),[],'all'));shading interp;
% title('Proposed');
% 
% figure(4);
% subplot(1,2,1);semilogy(reshape(stateP.NMSE,[],1));hold on;xlabel('Iteration',Interpreter='tex');
% subplot(1,2,2);semilogy(stateP.NMSE(end,:));hold on;xlabel('Iteration',Interpreter='tex');
% % switch inputParam.convBreaker
% %     case true
% %         semilogy(reshape(stateP.NMSE,[],1));hold on;xlabel('Iteration',Interpreter='tex');
% %     case false
% %         semilogy(stateP.NMSE(end,:));hold on;xlabel('Iteration',Interpreter='tex');
% % end
% % semilogy(reshape(nmse1(nmse1>0),[],1));hold on;xlabel('Iteration',Interpreter='latex');
% ylabel('NMSE',Interpreter='tex');
% title('Proposed');
% % semilogy(vecnorm((stateP.xUAVRec - xUAV),2,1));
% 
% % % figure(4);title('Sparsity Over Iterations',Interpreter='tex');
% % % semilogy(reshape(abs(spa),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % % semilogy(reshape(nmse2,[],1));legend('NMSE1','NMSE2',Location='best');
% % % % figure(5);
% % % % semilogy(reshape(abs(state.mu),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % % title('\mu',Interpreter='tex');
% % % figure(6);
% % % semilogy(reshape(abs(state.gammaP),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % title('\gamma_p',Interpreter='tex');
% % % figure(7);
% % % semilogy(reshape(abs(state.lambda),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % title('\lambda',Interpreter='tex');
% % % 
% % % figure(8);
% % % semilogy(reshape(abs(state.eta1),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % % title('\eta_1',Interpreter='tex');
% % % hold on;
% % % semilogy(reshape(abs(state.gamma1),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % % title('\eta_1',Interpreter='tex');
% % % semilogy(reshape(abs(state.gamma2),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % legend('\eta_1','\gamma_1','\gamma_2',Interpreter='tex');
% % % 
% % % figure(9);
% % % semilogy(reshape(abs(state.gammaOmega),[],1));hold on;xlabel('Iteration',Interpreter='tex');
% % % title('\gamma_\omega',Interpreter='tex');
%%
% figure(10);
% for ii = 1:3*niter
%     Itmp = vec2imag(state.xRec(:,ii));
% %     surf(abs(Itmp)/max(abs(Itmp),[],'all'));shading interp;
%     surf(abs(Itmp));shading interp;
%     strTitle = sprintf('%d EM iteration, %d VAMP iteration',...
%         37 + floor(ii/niter),mod(ii,niter)+1);
%     title(strTitle);
%     pause(0.5);
% end
% 
% figure(11);
% for ii = 1:inputParam.EMiter*inputParam.niter
%     plot(abs(state.d(:,ii)));shading interp;
%     strTitle = sprintf('%d EM iteration, %d VAMP iteration',...
%         ceil(ii/inputParam.niter),mod(ii,inputParam.niter)+1);
%     title(strTitle);
%     ylim([0,1]);
%     pause(0.3);
% end

% % figure(11);
% % for ii = 1:3*niter
% %     Itmp = vec2imag(state.mRec(:,ii));
% % %     surf(abs(Itmp)/max(abs(Itmp),[],'all'));shading interp;
% %     surf(abs(Itmp));shading interp;
% %     strTitle = sprintf('%d EM iteration, %d VAMP iteration',...
% %         37 + floor(ii/niter),mod(ii,niter)+1);
% %     title(strTitle);
% %     pause(0.5);
% % end
% figure(12);
% % for ii = 1:3*niter
% %     Itmp = vec2imag(abs(state.xRec(:,ii))-abs(state.mRec(:,ii)));
% % %     surf(abs(Itmp)/max(abs(Itmp),[],'all'));shading interp;
% %     surf((Itmp));shading interp;
% %     strTitle = sprintf('%d EM iteration, %d VAMP iteration',...
% %         37 + floor(ii/niter),mod(ii,niter)+1);
% %     title(strTitle);
% %     pause(0.3);
% % end
% % figure(13);
% % for ii = 1:inputParam.EMiter
% %     Itmp = vec2imag(state.muP(:,ii));
% % %     surf(abs(Itmp)/max(abs(Itmp),[],'all'));shading interp;
% %     surf(abs(Itmp));shading interp;
% %     strTitle = sprintf('%d EM iteration'...
% %         ,ii);
% %     title(strTitle);
% %     pause(0.5);
% % end

% % figure(14);
% % for ii = 1:3*niter
% %     Itmp = vec2imag(state.xrRec(:,ii));
% % %     surf(abs(Itmp)/max(arbs(Itmp),[],'all'));shading interp;
% %     surf(abs(Itmp));shading interp;
% %     switch mod(ii,3)
% %         case 0
% %             strTitle = sprintf('$r1$, %d-th iteration',ceil(ii/3));
% %         case 1
% %             strTitle = sprintf('$x1$, %d-th iteration',ceil(ii/3));
% %         case 2
% %             strTitle = sprintf('$r2$, %d-th iteration',ceil(ii/3));
% %     end
% %     title(strTitle,Interpreter="latex");
% %     pause(0.3);
% % end
%% Laplacian-VAMP
% xest3 = laplavamp(rEq,accuratePhiEq,1/sigma2/1000^2,niter);
% xest3 = xest3(1:length(xest3)/2)+1j*xest3(length(xest3)/2+1:end);
% I3 = vec2imag(xest3);
% 
% 
% % xest2 = amp(rEq,inaccuratePhiEq,1500,niter);
% % xest2 = xest2(1:length(xest2)/2)+1j*xest2(length(xest2)/2+1:end);
% % xest2 = amp(r,inaccuratePhi,niter);
% xest4 = laplavamp(rEq,inaccuratePhiEq,1/sigma2/1000^2,niter);
% xest4 = xest4(1:length(xest4)/2)+1j*xest4(length(xest4)/2+1:end);
% I4 = vec2imag(xest4);
% 
% figure(3);
% % subplot(1,2,1);
% surf(abs(I3)/max(abs(I3),[],'all'));shading interp;
% % subplot(1,2,2);
% figure(4);
% surf(abs(I4)/max(abs(I4),[],'all'));shading interp;



%% Orthogonal Match Pursuit
% xest5 = omp(r,accuratePhi,9);
% xest5 = xest5(Nt*Nr+1:end);
% 
% I5 = vec2imag(xest5);
% figure(5);
% surf(abs(I5)/max(abs(I5),[],'all'));shading interp;
% 
% xest6 = omp(r,inaccuratePhi,9);
% xest6 = xest6(Nt*Nr+1:end);
% I6 = vec2imag(xest6);
% figure(6);
% surf(abs(I6)/max(abs(I6),[],'all'));shading interp;

%% FISTA
% [xest7,err]=fista(r,accuratePhi,0.1,500);
% I7 = vec2imag(xest7);
% figure(7);
% surf(abs(I7)/max(abs(I7),[],'all'));shading interp;