%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Convergence Behavior\n');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle gg xRange dd
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)

deltaF = 120e+3;   % subcarrier spacing
Ns = 512;
BW = deltaF*Ns;  
T = 1 / deltaF;   % symbol duration
Tcp = T / 4;      % cyclic prefix duration
Ts = T + Tcp;     %   total symbol duration
sigmaPos = 1/sqrt(3);
sigma2 = 10^(-50/10);
xRange = 40;
niter = 20;
Na = 16;
Nv = 1;
Nsize = 21;
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
% xR = [-42; 0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = 1;
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
pilotSubcarrInd = 0:32:Ns-1;
dataSubcarrInd = [];
scaleFac = 5;
for ii = 1:scaleFac
dataSubcarrInd = [dataSubcarrInd , pilotSubcarrInd + ii];
end
dataSubcarrInd = sort(dataSubcarrInd);
dataNum = length(dataSubcarrInd);

% dataNum = length(dataSubcarrInd);

% pilot = ones(length(pilotSubcarrInd),1);
% symbols = ones(dataNum,1);
% constell = 1/sqrt(2)*[1+1j,-1+1j,-1-1j,1-1j];
constell(1).Cons = 1/sqrt(2)*[1+1j,-1+1j,-1-1j,1-1j];
constell(2).Cons = [1/sqrt(2)*[1+1j,-1+1j,...
                             -1-1j,1-1j],1+0j,0+1j,-1+0j,0-1j];
constell(3).Cons = 1/sqrt(10)*[-3+3j,-1+3j,1+3j,3+3j,...
                       -3+1j,-1+1j,1+1j,3+1j,...
                       -3-1j,-1-1j,1-1j,3-1j,...
                       -3-3j,-1-3j,1-3j,3-3j];
% xUAV = [50,0,0;0,50,-50;60,60,60];
Nt = 6;
Nsim = 10;
% xUAV = [50,30,0,0,0,0;0,0,30,-30,50,-50;60,60,60,60,60,60];
% xUAV = 100*rand(2,Nt)-50;
% xUAV = [xUAV;ones(1,Nt)*60];
% xUAV = [50;0;60];


% Nt = size(xUAV,2);
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
% xUAVEst = xUAV + sigmaPos*(rand(size(xUAV))-1/2);
% xUAVEst = xUAV;

% vecTO = 0*1/BW*(rand(Nr,1) - 1/2)/10;
vecTO = 0/BW*(rand-1/2)*ones(Nr,1)/2; % offset is ignored
% vecTO = 1/BW*(rand-1/2);

dArray = c0/fc/2;
biasAngle = [0]';                                                         %degree
% biasAngle = [90]';                                                      %degree

% biasAngle = zeros(Nr,2); 
% vecTO = zeros(Nr,1);
vecTOEst = zeros(size(vecTO));
% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
xTarget = zeros(3,16);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation
% Phip = zeros(Na*Nv*length(pilotSubcarrInd),Nsize^2+1,Nt);
% Phid = zeros(Na*Nv*length(dataSubcarrInd),Nsize^2+1,Nt);
% yp = zeros(size(Phip,1),Nt);
% yd = zeros(size(Phid,1),Nt);
% 
% numTarget = 13;
% [targetInd] = randperm(Nsize^2,numTarget);
% PhidwithSym = Phid;
% d = zeros(dataNum,Nt);
% x = zeros(Nsize^2+1,Nt);
eyeS = eye(Na*Nv);
% for tt = 1:Nt
%     Phip(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,pilotSubcarrInd,pilot,  vecTO);
%     Phid(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,symbols,vecTO);
% 
%     bits = randi([0,1],dataNum,2);
%     d(:,tt) = (2*bits(:,1)-1) + 1j*(2*bits(:,2)-1);
%     
% %     PhidwithSym(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,d(:,tt),vecTO);
%     PhidwithSym(:,:,tt) = repmat(kron(ones(Na*Nv,1),d(:,tt)),1,size(Phid,2)).*Phid(:,:,tt);
%     if tt > 1
%         x(1+targetInd,tt) = 0.7*x(1+targetInd,tt-1) + 0.3/sqrt(2)*(randn(numTarget,1) + 1j*randn(numTarget,1));
%     else
%         x(1+targetInd,tt) = 1 + 1/sqrt(2)*(randn(numTarget,1) + 1j*randn(numTarget,1));
%     end
%     x(1,tt) = 5 + 10/sqrt(2)*(randn + 1j*randn);
%     yp(:,tt) = Phip(:,:,tt)*x(:,tt) + (normrnd(0,sqrt(sigma2/2),size(Phip,1),1) +...
%                     1j*normrnd(0,sqrt(sigma2/2),size(Phip,1),1));
%     yd(:,tt) = PhidwithSym(:,:,tt)*x(:,tt) + (normrnd(0,sqrt(sigma2/2),size(Phid,1),1) +...
%                     1j*normrnd(0,sqrt(sigma2/2),size(Phid,1),1));
% end
% gg = x;
UAVarray = 1:2:16;
load('dataset\FusNum\Phid.mat');
load('dataset\FusNum\Phip.mat')
load('dataset\FusNum\PhidGenie.mat');
load('dataset\FusNum\PhipGenie.mat');

load('dataset\FusNum\xUAV.mat');

load("dataset\FusNum\d.mat");
load("dataset\FusNum\noisedatap.mat");
load("dataset\FusNum\noisedatad.mat");
load("dataset\FusNum\xOG.mat");
load("dataset\FusNum\zp.mat");
load("dataset\FusNum\zd.mat");
% load("dataset\FusNum\derivD.mat");
% load("dataset\FusNum\derivP.mat");
load("dataset\FusNum\sMGParam.mat");
load('dataset\FusNum\grid.mat');
% sMGParam.derivD = derivD;
% sMGParam.derivP = derivP;

% d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);

gg = x(:,UAVarray);
snr = 5; nMod = 3;
dd = d(:,UAVarray);

% figCount = 1;
signalPower = mean(vecnorm([zp;zd(:,:)],2).^2)/(size(zp,1) + size(zd,1));
signalPowerp = mean(vecnorm(zp,2).^2)/size(zp,1); 
signalPowerd = mean(vecnorm(zd(:,:),2).^2)/size(zd,1);

% snr = 10; nMod = 3;
for nsim = 1:5
    yp = zp + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,nsim);
    yd = zd(:,:) + sqrt(signalPowerd/(10^(snr/10)))*noisedatad(:,:,nsim);
    %% VBI
    % clear stateSBL
    figCount = 1;
    inputParam.pilotSubcarrInd = pilotSubcarrInd;
    inputParam.dataSubcarrInd  = dataSubcarrInd;
    inputParam.constell = constell(nMod).Cons;
    inputParam.gammaOmega = 1e5;
    inputParam.dataNum = dataNum;
    inputParam.eyeS = eyeS;
    inputParam.SV = 0;  
    inputParam.Niter = 80;
    inputParam.NiterA = 20;
    inputParam.JESD = 1;
    inputParam.xUAV = xUAV(:,UAVarray);
    % inputParam.xR = xR;
    % inputParam.NiterA = 20;
    inputParam.lambda0 = 0.90;
    inputParam.lambdaS = 0.005;
    % inputParam.lambda0 = 0.8;
    inputParam.convBreak = 0;
    
    [xestVBI,destVBI,stateVBI] = VBI2(yp(:,UAVarray),yd(:,UAVarray),Phip(:,:,UAVarray),Phid(:,:,UAVarray),inputParam);
    nmseVBI = reshape(stateVBI.NMSE,[],1);
    nmseVBI = nmseVBI(nmseVBI~=0);
    % NMSEVBI = vec(stateVBI.NMSE(end,:));
%     NMSEVBI = find4Conv(stateVBI.NMSE); 
    % NMSEVBI = [nmseVBI(1),NMSEVBI];
    dDecodedVBI = mleDecoder(destVBI,constell(nMod).Cons);
    SERVBI = sum(dDecodedVBI~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    NMSEVBI(nsim,:) = [nmseVBI(1),find4Conv(stateVBI.NMSE)];
    serVBI(nsim,:) = stateVBI.SER;
    
    
    %% VBIWOD
    % clear stateSBL
    % figCount = 1;
    inputParam.pilotSubcarrInd = pilotSubcarrInd;
    inputParam.dataSubcarrInd  = dataSubcarrInd;
    inputParam.constell = constell(nMod).Cons;
%     inputParam.gammaOmega = 1e6;
    inputParam.dataNum = dataNum;
    inputParam.eyeS = eyeS;
    inputParam.SV = 0;  
    inputParam.Niter = 80;
    inputParam.NiterA = 20;
    inputParam.JESD = 0    ;
    inputParam.xUAV = xUAV(:,UAVarray);
    % inputParam.xR = xR;
    % inputParam.NiterA = 20;
    inputParam.lambda0 = 0.90;
    inputParam.lambdaS = 0.005;
    % inputParam.lambda0 = 0.8;
    inputParam.convBreak = 0;
    
    [xestVBIWOD,~,stateVBIWOD] = VBI2(yp(:,UAVarray),yd(:,UAVarray),Phip(:,:,UAVarray),Phid(:,:,UAVarray),inputParam);
    nmseVBIWOD = reshape(stateVBIWOD.NMSE,[],1);
    nmseVBIWOD = nmseVBIWOD(nmseVBIWOD~=0);
    % NMSEVBI = vec(stateVBI.NMSE(end,:));
%     NMSEVBIWOD = find4Conv(stateVBIWOD.NMSE); 
    NMSEVBIWOD(nsim,:) = [nmseVBIWOD(1),find4Conv(stateVBIWOD.NMSE)];
    %% Baseline1 (SVSBL)
    % figCount = 1;
    % inputParam.constell = constell;
    inputParam.pilotSubcarrInd = pilotSubcarrInd;
    inputParam.dataSubcarrInd  = dataSubcarrInd;
    % inputParam.gammaOmega = 1e5;
    inputParam.dataNum = dataNum;
    inputParam.eyeS = eyeS;
    % inputParam.Niter = 20;
    % inputParam.NiterA = 20;
    inputParam.JESD = 1;
    inputParam.SV = 1;
    inputParam.sMGParam = sMGParam;
    % inputParam.xR = xR;
    inputParam.armijoSigma = 0.000;
    inputParam.armijoBeta = 0.4;
    inputParam.attenFac = 0.96;

    % inputParam.NiterA = 20;
    % inputParam.lambda0 = 0.5;
    % inputParam.lambdaS = 0.002;
    
    % inputParam.lambda0 = 0.8;
    % inputParam.convBreak = 0;
    % for tt = 1:Nt
    %     gg = x(:,tt);
    %     [xestSV(:,tt),destSV(:,tt),stateSV(tt)] = SVVBI(yp(:,tt),yd(:,tt),Phip(:,:,tt),Phid(:,:,tt),inputParam);
    % end
    tic
    [xestSV,destSV,stateSV] = SVVBIOG(yp(:,UAVarray),yd(:,UAVarray),Phip(:,:,UAVarray),Phid(:,:,UAVarray),inputParam);
    toc
    nmseSV = reshape(stateSV.NMSE,[],1);
    nmseSV = nmseSV(nmseSV~=0);
    % NMSESV = vec(stateSV.NMSE(end,:));
    NMSESV(nsim,:) = find4Conv(stateSV.NMSE);
    serSV(nsim,:) = stateSV.SER;
    
    % figure(figCount);figCount = figCount + 1;
    % for ii =1:Nt
    %     subplot(1,Nt,ii);
    %     stem(abs(xestSV(:,ii)),'o');
    %     hold on;stem(abs(x(:,ii)),'x');
    % end
    % nmseSV = zeros(size((stateSV(1).NMSE)));
    % for tt = 1:Nt
    %     nmseSV = nmseSV + (stateSV(tt).NMSE);
    % end
    % nmseSV = nmseSV/Nt;
    % % figure(figCount);figCount = figCount + 1;
    % NMSESV = vec(nmseSV(end,:)); 
    % nmseSV = reshape(nmseSV,[],1);
    
    % semilogy(nmseSV);
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     semilogy(vec(stateSV(tt).NMSE));hold on;
    % end
    if inputParam.JESD
        dDecodedSV = mleDecoder(destSV,constell(nMod).Cons);
        SERSV  =  sum(dDecodedSV~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    end
       
    % 
    % figure(figCount);figCount = figCount + 1;
    % nmseAL = reshape(stateAL(1).NMSE,[],1);
    
    % semilogy(nmseAL);
    
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     subplot(1,Nt,tt)
    %     semilogy(stateAL(tt).rho);
    % end
    %% Baseline2
    % gg = x;
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.sMGParam = sMGParam;
    % inputParam.decoderFlag = 0;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR;
    % inputParam.JESD = 1;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % [xestAL,destAL,stateAL] = ALESDOG(yp,yd,Phip,Phid,inputParam);
    % toc
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:Nt
    % %     subplot(1,Nt,ii);
    % %     stem(abs(xestAL(:,ii)),'o');
    % %     hold on;stem(abs(x(:,ii)),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % nmseAL = reshape(stateAL.NMSE,[],1);
    % nmseAL = nmseAL(nmseAL~=0);
    % % NMSEAL = vec(stateAL.NMSE(end,:));
    % NMSEAL = find4Conv(stateAL.NMSE);
    % % semilogy(nmseAL);
    % 
    % SERAL  =  sum(destAL~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    
    %% Baseline 3
%     gg = x;
    inputParam.pilotSubcarrInd = pilotSubcarrInd;
    inputParam.dataSubcarrInd  = dataSubcarrInd;
    % inputParam.gammaOmega = 1./sigma2;
    inputParam.dataNum = dataNum;
    inputParam.sMGParam = sMGParam;
    inputParam.eyeS = eyeS;
    % inputParam.Niter = 10;
    inputParam.xUAV = xUAV(:,UAVarray);
    % inputParam.xR = xR; 
    inputParam.JESD = 0;
    inputParam.armijoSigma = 0.000;
    inputParam.armijoBeta = 0.4;
    inputParam.lambda0 = 0.90;
    inputParam.lambdaS = 0.005;
    % inputParam.NiterA = 10;
    % inputParam.lambda0 = 0.05;
    inputParam.convBreak = 0; 
    tic
    [xestWOD,destWOD,stateWOD] = MVVBIOG(yp(:,UAVarray),yd(:,UAVarray),Phip(:,:,UAVarray),Phid(:,:,UAVarray),inputParam);
    % for tt = 1:Nt
    %     Phidd(:,:,tt) = repmat(kron(ones(Na*Nv,1),d(:,tt,nMod)),1,size(Phid,2)).*PhidGenie(:,:,tt);
    % end 
    % [xestWOD11,destWOD11,stateWOD11] = MVVBI([yp;yd],yd,[PhipGenie;Phidd],Phid,inputParam);
    % [xestWOD22,destWOD22,stateWOD22] = MVVBI(yd,yd,Phidd,Phid,inputParam);
    % [xestWOD33,destWOD33,stateWOD33] = MVVBI(yp,yd,PhipGenie,Phid,inputParam);  
    toc
    
    % figure(figCount);figCount = figCount + 1;
    % for ii =1:1
    %     subplot(1,1,ii);
    %     stem(abs(xestWOD(:,ii)),'o');
    %     hold on;stem(abs(x(:,ii)),'x');
    % end
    % figure(figCount);figCount = figCount + 1;
    nmseWOD = reshape(stateWOD.NMSE,[],1);
    nmseWOD = nmseWOD(nmseWOD~=0);
    % NMSEWOD = vec(stateWOD.NMSE(end,:));
    NMSEWOD(nsim,:) = find4Conv(stateWOD.NMSE);
    % semilogy(nmseWOD);
    
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     subplot(1,Nt,tt);
    %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % end
    
    % figure(figCount);figCount = figCount + 1;
    % semilogy(nmseWOD);
    % for tt = 1:Nt
    %     stem(state.st(:,tt));hold on;
    % end
    % set(gca,'YScale','log');
    
    % vecnorm(xest-gg,2).^2
    
    
    % %% Baseline 3（OG）
    % gg = x;
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.sMGParam = sMGParam;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR; 
    % inputParam.JESD = 0;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % [xestWOD2,destWOD2,stateWOD2] = MVVBIOG2(yp,yd,Phip,Phid,inputParam);
    % % [xestWOD,destWOD,stateWOD] = MVVBI(yp,yd,Phip,Phid,inputParam);
    % toc
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:1
    % %     subplot(1,1,ii);
    % %     stem(abs(xestWOD(:,ii)),'o');
    % %     hold on;stem(abs(x(:,ii)),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % nmseWOD2 = reshape(stateWOD2.NMSE,[],1);
    % nmseWOD2 = nmseWOD2(nmseWOD2~=0);
    % % NMSEWOD = vec(stateWOD.NMSE(end,:));
    % NMSEWOD2 = find4Conv(stateWOD2.NMSE);
    % semilogy(nmseWOD);
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     subplot(1,Nt,tt);
    %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % end
    
    % figure(figCount);figCount = figCount + 1;
    % semilogy(nmseWOD);
    % for tt = 1:Nt
    %     stem(state.st(:,tt));hold on;
    % end
    % set(gca,'YScale','log');
    
    % vecnorm(xest-gg,2).^2
    % %% Baseline 3（OG3）
    % gg = x;
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.sMGParam = sMGParam;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR; 
    % inputParam.JESD = 0;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % [xestWOD3,destWOD3,stateWOD3] = MVVBIOG3(yp,yd,Phip,Phid,inputParam);
    % % [xestWOD,destWOD,stateWOD] = MVVBI(yp,yd,Phip,Phid,inputParam);
    % toc
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:1
    % %     subplot(1,1,ii);
    % %     stem(abs(xestWOD(:,ii)),'o');
    % %     hold on;stem(abs(x(:,ii)),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % nmseWOD3 = reshape(stateWOD3.NMSE,[],1);
    % nmseWOD3 = nmseWOD3(nmseWOD3~=0);
    % % NMSEWOD = vec(stateWOD.NMSE(end,:));
    % NMSEWOD3 = find4Conv(stateWOD3.NMSE);
    % % semilogy(nmseWOD);
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     subplot(1,Nt,tt);
    % %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % % end
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % semilogy(nmseWOD);
    % % for tt = 1:Nt
    % %     stem(state.st(:,tt));hold on;
    % % end  
    % % set(gca,'YScale','log');
    % 
    % % vecnorm(xest-gg,2).^2
    % %% Proposed Algorithm(OG3)
    % dbstop if error;
    % gg = x;
    % dd = d(:,UAVarray);
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.sMGParam = sMGParam;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR;
    % inputParam.JESD = 1;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % [xestOG3,destOG3,stateOG3] = MVVBIOG3(yp,yd,Phip,Phid,inputParam);
    % toc
    % % [xest,dest,state] = MVVBIOG(yp,yd,Phip,Phid,inputParam);
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:Nt
    % %     subplot(1,Nt,ii);
    % %     stem(abs(xestOG(:,ii)),'o'); 
    % %     hold on;stem(abs(x(:,ii) ),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % % nmse = reshape(state.NMSE,[],1);
    % % nmse = nmse(nmse~=0);
    % nmseOG3 = reshape(stateOG3.NMSE,[],1);
    % nmseOG3 = nmseOG3(nmseOG3~=0);
    % % NMSE = vec(state.NMSE(end,:));
    % NMSEOG3 = find4Conv(stateOG3.  MSE);
    % % semilogy(nmse);
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     subplot(1,Nt,tt);
    % %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % % end
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     stem(state.st(:,tt));hold on;
    % % end
    % % set(gca,'YScale','log');
    % 
    % % vecnorm(xest-gg,2).^2
    % if inputParam.JESD
    %     dDecodedOG3 = mleDecoder(destOG3,constell(nMod).Cons);
    %     SERProposedOG3  =  sum(dDecodedOG3~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    % end
    %% Proposed Algorithm(OG)
    dbstop if error;
    % gg = x(:,1);
    % dd = d(:,1,nMod);
%     gg = x;
%     dd = d(:,UAVarray);
    inputParam.pilotSubcarrInd = pilotSubcarrInd;
    inputParam.dataSubcarrInd  = dataSubcarrInd;
    % inputParam.gammaOmega = 1./sigma2;
    inputParam.dataNum = dataNum;
    inputParam.eyeS = eyeS;
    inputParam.Niter = 80;
    inputParam.sMGParam = sMGParam;
    inputParam.xUAV = xUAV(:,UAVarray);
    % inputParam.xR = xR; 
    inputParam.JESD = 1;
    inputParam.armijoSigma = 0.01;
    inputParam.armijoBeta = 0.4;
    inputParam.lambda0 = 0.90;
    inputParam.lambdaS = 0.005;
    % inputParam.NiterA = 10; 
    % inputParam.lambda0 = 0.05;
    inputParam.convBreak = 0;  
    tic
    [xestOG,destOG,stateOG] = MVVBIOG(yp(:,UAVarray),yd(:,UAVarray),Phip(:,:,UAVarray),Phid(:,:,UAVarray),inputParam);
%     [xestOG,destOG,stateOG] = MVVBI(yp,yd,PhipGenie,PhidGenie,inputParam);
    
    toc  
    % [xest,dest,state] = MVVBIOG(yp,yd,Phip,Phid,inputParam);    
    
    % figure(figCount);figCount = figCount + 1;
    % for ii =1:Nt
    %     subplot(1,Nt,ii);
    %     stem(abs(xestOG(:,ii)),'o'); 
    %     hold on;stem(abs(x(:,ii) ),'x');
    % end
    % figure(figCount);figCount = figCount + 1;
    % nmse = reshape(state.NMSE,[],1);
    % nmse = nmse(nmse~=0);
    nmseOG = reshape(stateOG.NMSE,[],1);
    nmseOG = nmseOG(nmseOG~=0);
    % NMSE = vec(state.NMSE(end,:));
    NMSEOG(nsim,:) = find4Conv(stateOG.NMSE);
    serOG(nsim,:) = stateOG.SER;
    % semilogy(nmse);
    
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     subplot(1,Nt,tt);
    %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % end
    
    % figure(figCount);figCount = figCount + 1;
    % for tt = 1:Nt
    %     stem(state.st(:,tt));hold on;
    % end
    % set(gca,'YScale','log');
    
    % vecnorm(xest-gg,2).^2
    if inputParam.JESD
        dDecodedOG = mleDecoder(destOG,constell(nMod).Cons);
        SERProposedOG  =  sum(dDecodedOG~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    end
    % if inputParam.JESD
    %     dDecodedOG = mleDecoder(destOG,constell(nMod).Cons);
    %     SERProposedOG  =  sum(dDecodedOG~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    % end
      
    %% Proposed Algorithm(OG2)
    % dbstop if error;
    % gg = x;
    % dd = d(:,UAVarray);
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.sMGParam = sMGParam;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR;
    % inputParam.JESD = 1;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % [xestOG2,destOG2,stateOG2] = MVVBIOG2(yp,yd,Phip,Phid,inputParam);
    % toc
    % % [xest,dest,state] = MVVBIOG(yp,yd,Phip,Phid,inputParam);
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:Nt
    % %     subplot(1,Nt,ii);
    % %     stem(abs(xestOG(:,ii)),'o'); 
    % %     hold on;stem(abs(x(:,ii) ),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % % nmse = reshape(state.NMSE,[],1);
    % % nmse = nmse(nmse~=0);
    % nmseOG2 = reshape(stateOG2.NMSE,[],1);
    % nmseOG2 = nmseOG2(nmseOG2~=0);
    % % NMSE = vec(state.NMSE(end,:));
    % NMSEOG2 = find4Conv(stateOG2.NMSE);
    % % semilogy(nmse);
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     subplot(1,Nt,tt);
    % %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % % end
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     stem(state.st(:,tt));hold on;
    % % end
    % % set(gca,'YScale','log');
    % 
    % % vecnorm(xest-gg,2).^2
    % if inputParam.JESD
    %     dDecodedOG2 = mleDecoder(destOG2,constell(nMod).Cons);
    %     SERProposedOG2  =  sum(dDecodedOG2~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    % end
    % % if inputParam.JESD
    % %     dDecodedOG = mleDecoder(destOG,constell(nMod).Cons);
    % %     SERProposedOG  =  sum(dDecodedOG~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    % % end
    %%
    
    
    %% Proposed Algorithm
    % dbstop if error;
    % gg = x;
    % dd = d(:,UAVarray);
    % 
    % inputParam.pilotSubcarrInd = pilotSubcarrInd;
    % inputParam.dataSubcarrInd  = dataSubcarrInd;
    % % inputParam.gammaOmega = 1./sigma2;
    % inputParam.dataNum = dataNum;
    % inputParam.eyeS = eyeS;
    % % inputParam.Niter = 10;
    % inputParam.xUAV = xUAV;
    % inputParam.xR = xR;
    % inputParam.JESD = 1;
    % % inputParam.NiterA = 10;
    % % inputParam.lambda0 = 0.05;
    % % inputParam.convBreak = 1;
    % tic
    % % [xest,dest,state] = MVVBI(yp,yd,Phip,Phid,inputParam);
    % [xest,dest,state] = MVVBI(yp,yd,PhipGenie,PhidGenie,inputParam);
    % 
    % toc
    % % figure(figCount);figCount = figCount + 1;
    % % for ii =1:Nt
    % %     subplot(1,Nt,ii);
    % %     stem(abs(xest(:,ii)),'o');
    % %     hold on;stem(abs(x(:,ii)),'x');
    % % end
    % % figure(figCount);figCount = figCount + 1;
    % % nmse = reshape(state.NMSE,[],1);
    % % nmse = nmse(nmse~=0);
    % nmse = reshape(state.NMSE,[],1);
    % nmse = nmse(nmse~=0);
    % % NMSE = vec(state.NMSE(end,:));
    % NMSE = find4Conv(state.NMSE);
    % % semilogy(nmse);
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     subplot(1,Nt,tt);
    % %     semilogy(state.rho(:,tt));hold on;semilogy(stateSV(tt).rho);
    % % end
    % 
    % % figure(figCount);figCount = figCount + 1;
    % % for tt = 1:Nt
    % %     stem(state.st(:,tt));hold on;
    % % end
    % % set(gca,'YScale','log');
    % 
    % % vecnorm(xest-gg,2).^2
    % if inputParam.JESD
    %     dDecoded = mleDecoder(dest,constell(nMod).Cons);
    %     SERProposed  =  sum(dDecoded~=d(:,UAVarray),'all')/numel(d(:,UAVarray))
    % end
end
%% Plot
modStr = strings(3,1);
modStr(1) = '4QAM';modStr(2) = '8PSK';modStr(3) = '16QAM';

figure(figCount);figCount = figCount + 1;

% semilogy(nmseOG,'-o');hold on;semilogy(nmseVBI,':>');semilogy(nmseSV,':s');semilogy(nmseVBIWOD,'--<');semilogy(nmseWOD,'-.d');xlim([1 100]);
semilogy(nmseOG,'-o');hold on;semilogy(nmseVBI,':>');semilogy(nmseSV,':s');
% semilogy(nmseAL,':p');
semilogy(nmseVBIWOD,'--<');semilogy(nmseWOD,'-.d');xlim([1 100]);
% legend('OGMVSBL','VBI','SVVBL','ALESD','SBL','MVSBLWOD',Interpreter='tex',FontName='Times New Roman',FontSize=10);

legend('MVSBL','VBI','SVVBL','SBL','MVSBLWOD',Interpreter='tex',FontName='Times New Roman',FontSize=10);
ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
grid on;
% nameStr = strcat('fig\figureConvergenceFullOG',modStr(nMod),'.fig');
% savefig(nameStr);
%%
Niter = 80;
figure(figCount);figCount = figCount + 1;
downSampFac = 3;

% semilogy(NMSE,'-o');
% hold on;
semilogy(1:downSampFac:Niter,mean(serVBI(:,1:downSampFac:Niter),1),'--<',LineWidth=1.5,Color = [255,0,0]/255);hold on;
semilogy(1:downSampFac:Niter,mean(serSV(:,1:downSampFac:Niter),1),':>',LineWidth=1.5,Color = [255,128,0]/255);
semilogy(1:downSampFac:Niter,mean(serOG(:,1:downSampFac:Niter),1),'-o',LineWidth=1.5,Color = [0,0,255]/255);
% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend('VBI','SVVBI','MVVBI',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);
ylim([0.2e-2,1]);
ylabel('SER',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
grid on;
nameStr = strcat('fig\figureSERConvergenceOG',modStr(nMod),'.fig');
savefig(nameStr);

%%
figure(figCount);figCount = figCount + 1;
downSampFac = 3;

% semilogy(NMSE,'-o');
% hold on;
semilogy(1:downSampFac:Niter,mean(NMSEVBIWOD(:,1:downSampFac:Niter),1),':s',LineWidth=1.5,Color = [255,0,255]/255);hold on;
semilogy(1:downSampFac:Niter,mean(NMSEVBI(:,1:downSampFac:Niter),1),'--<',LineWidth=1.5,Color = [255,0,0]/255);
semilogy(1:downSampFac:Niter,mean(NMSESV(:,1:downSampFac:Niter),1),':>',LineWidth=1.5,Color = [255,128,0]/255);
% semilogy(1:downSampFac:Niter,NMSEAL(1:downSampFac:Niter),':d');
semilogy(1:downSampFac:Niter,mean(NMSEWOD(:,1:downSampFac:Niter),1),'-.d',LineWidth=1.5,Color = [0,255,0]/255);
semilogy(1:downSampFac:Niter,mean(NMSEOG(:,1:downSampFac:Niter),1),'-o',LineWidth=1.5,Color = [0,0,255]/255);
% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend('SBL','VBI','SVVBI','MVVBIWOD','MVVBI',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);
ylim([0.2e-2,1]);
ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
grid on;
nameStr = strcat('fig\figureConvergenceOG2',modStr(nMod),'.fig');
savefig(nameStr);

%%
str = datestr(now);
str = strrep(str,' ','');
str = strrep(str,'-','');
str = strrep(str,':','');

str = strcat('results\',str,'mainOGConvtest',modStr(nMod),'.mat');
save(str,'xestWOD','xestOG','xestSV','xestVBI','xestVBIWOD',...
    'destSV','destOG','destVBI','NMSEVBI','NMSEVBIWOD','NMSEOG',...
    'NMSEWOD','NMSESV','nmseVBI','nmseVBIWOD','nmseOG','nmseWOD','nmseSV');
    save(str,'xestWOD','xestOG','xestSV','xestVBI','xestVBIWOD','xestAL',...
    'destSV','destOG','destVBI','destAL','NMSEVBI','NMSEVBIWOD','NMSEOG','NMSEWOD','NMSESV','NMSEAL','nmseVBI','nmseVBIWOD','nmseOG','nmseWOD','nmseSV','nmseAL');
fprintf('\n Convergence Behavior Simulation Finished!\n');
%