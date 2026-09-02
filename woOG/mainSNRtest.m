%% mainSNRtest
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Performance vs SNR\n');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle gg xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)

deltaF = 120e+3;   % subcarrier spacing
Ns = 256;
BW = deltaF*Ns;  
T = 1 / deltaF;   % symbol duration
Tcp = T / 4;      % cyclic prefix duration
Ts = T + Tcp;     % total symbol duration
sigmaPos = 1/sqrt(3);
sigma2 = 10^(-60/10);
xRange = 40;
niter = 20;
Na = 16;
Nv = 1;
Nsize = 21;
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
xR = [-22; -0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = size(xR,2);
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
pilotSubcarrInd = 0:15:Ns-1;
dataSubcarrInd = 0:1:Ns-1;
dataSubcarrInd = setdiff(dataSubcarrInd,pilotSubcarrInd);
% dataSubcarrInd = pilotSubcarrInd;
% dataNum = length(dataSubcarrInd);

dataNum = length(dataSubcarrInd);

pilot = ones(length(pilotSubcarrInd),1);
symbols = ones(dataNum,1);
modNum = 3;
constell(1).Cons = 1/sqrt(2)*[1+1j,-1+1j,-1-1j,1-1j];
constell(2).Cons = 1/sqrt(2)*[sqrt(2)+0j,1+1j,0+1j*sqrt(2),-1+1j,...
                             -sqrt(2)+0j,-1-1j,0-1j*sqrt(2),1-1j];
constell(3).Cons = 1/sqrt(10)*[-3+3j,-1+3j,1+3j,3+3j,...
                       -3+1j,-1+1j,1+1j,3+1j,...
                       -3-1j,-1-1j,1-1j,3-1j,...
                       -3-3j,-1-3j,1-3j,3-3j];
% xUAV = [50,0,0;0,50,-50;60,60,60];
Nt = 10;
% xUAV = [50,30,0,0,0,0;0,0,30,-30,50,-50;60,60,60,60,60,60];
% xUAV = 100*rand(2,Nt)-50;
% xUAV = [xUAV;ones(1,Nt)*60];
% xUAV = [50;0;60];

Nsim = 10;

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
SNR = [-10,-5,0,5,10,15];
snr = 10.^(SNR./10);

% SNR = [-20,-15,-10,-5,0];

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
load('dataset\1\Phid.mat');
load('dataset\1\Phip.mat');
load('dataset\1\xUAV.mat');

load("dataset\1\d.mat");
load("dataset\1\noisedatap.mat");
load("dataset\1\noisedatad.mat");
load("dataset\x.mat");
load("dataset\1\zp.mat");
load("dataset\1\zd.mat");
gg = x;
d(:,:,1) = mleDecoder(d(:,:,1),constell(1).Cons);
d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);
d(:,:,3) = mleDecoder(d(:,:,3),constell(3).Cons);

signalPower = mean(vecnorm([zp;zd(:,:,1)]).^2)/(size(zp,1) + size(zd,1));
signalPowerp = mean(vecnorm(zp).^2)/size(zp,1);
signalPowerd = mean(vecnorm(zd(:,:,1)).^2)/size(zd,1);

% signalPowerp = repmat(vecnorm(zp).^2/size(zp,1),size(zp,1),1);
% signalPowerd = repmat(vecnorm(zd).^2/size(zd,1),size(zd,1),1);

% figCount = 1;
datestr(now)
figCount = 1;
tic
%% SBL
clear stateSBL

inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;

% inputParam.gammaOmega = 1e6;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
inputParam.SV = 0;
inputParam.Niter = 10;
inputParam.NiterA = 20;
inputParam.JESD = 1;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
% inputParam.NiterA = 20;
inputParam.lambda0 = 0.5;
inputParam.lambdaS = 0.001;
% inputParam.lambda0 = 0.8;
inputParam.convBreak = 1;
for nMod = 1:modNum
    inputParam.constell = constell(nMod).Cons;
    for snrIdx =1:length(SNR)
        inputParam.gammaOmega = 1e6*(mean(signalPower)/snr(snrIdx))^(-1);
        for nsim = 1:Nsim
            yp = zp + sqrt(signalPowerp/snr(snrIdx)).*noisedatap(:,:,nsim);
            yd = zd(:,:,nMod) + sqrt(signalPowerd/snr(snrIdx)).*noisedatad(:,:,nsim);

            [xestVBI(:,:,snrIdx,nsim,nMod),destVBI(:,:,snrIdx,nsim,nMod),stateVBI(snrIdx,nMod)] = VBI2(yp,yd,Phip,Phid,inputParam);
            NMSEVBI(snrIdx,nsim,nMod) = mean(vecnorm(xestVBI(:,:,snrIdx,nsim,nMod) - x,2,1).^2./vecnorm(x,2,1).^2);
            dDecodedVBI(:,:,snrIdx,nsim,nMod) = mleDecoder(destVBI(:,:,snrIdx,nsim,nMod),constell(nMod).Cons);
            SERVBI(snrIdx,nsim,nMod) = sum(dDecodedVBI(:,:,snrIdx,nsim,nMod)~=d(:,:,nMod),'all')/numel(d(:,:,nMod));
            fprintf('VBI, SNR = %d dB, nsim = %d, Modulation Mode %d.\n',SNR(snrIdx),nsim,nMod);
        end
    end
end
nmseVBI = reshape(mean(NMSEVBI,2),[],modNum);
%% Baseline1 (SVSBL)

clear stateSV
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.constell = constell;
% inputParam.gammaOmega = 1e5;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
inputParam.SV = 1;
inputParam.Niter = 10;
inputParam.NiterA = 20;
inputParam.JESD = 1;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
% inputParam.NiterA = 20;
% inputParam.lambda0 = 0.5;
% inputParam.lambdaS = 0.002;
% inputParam.lambda0 = 0.8;
inputParam.convBreak = 1;
for nMod = 1:modNum
    inputParam.constell = constell(nMod).Cons;
    for snrIdx =1:length(SNR)
        inputParam.gammaOmega = 1e6*(mean(signalPower)/snr(snrIdx))^(-1);    
        for nsim = 1:Nsim
            yp = zp + sqrt(signalPowerp/snr(snrIdx)).*noisedatap(:,:,nsim);
            yd = zd(:,:,nMod) + sqrt(signalPowerd/snr(snrIdx)).*noisedatad(:,:,nsim);
    %         for tt = 1:Nt
    %             inputParam.gammaOmega = ((size(Phid,1)+size(Phip,1))*snr(snrIdx))./(norm([yp(:,tt);yd(:,tt)])^2);
    %             gg = zeros(size(gg));
    %             [xestSV(:,tt,snrIdx,nsim),destSV(:,tt,snrIdx,nsim),~] = SVVBI(yp(:,tt),yd(:,tt),Phip(:,:,tt),Phid(:,:,tt),inputParam);
    %         end
    %         inputParam.gammaOmega = ((size(Phid,1)+size(Phip,1))*snr(snrIdx))./(norm([yp(:,tt);yd(:,tt)])^2);
            [xestSV(:,:,snrIdx,nsim,nMod),destSV(:,:,snrIdx,nsim,nMod),stateSV(snrIdx,nMod)] = SVVBI(yp,yd,Phip,Phid,inputParam);
            NMSESV(snrIdx,nsim,nMod) = mean(vecnorm(xestSV(:,:,snrIdx,nsim,nMod) - x,2,1).^2./vecnorm(x,2,1).^2);
            dDecodedSV(:,:,snrIdx,nsim,nMod) = mleDecoder(destSV(:,:,snrIdx,nsim,nMod),constell(nMod).Cons);
            SERSV(snrIdx,nsim,nMod) = sum(dDecodedSV(:,:,snrIdx,nsim,nMod)~=d(:,:,nMod),'all')/numel(d(:,:,nMod));
            fprintf('Baseline 1, SNR = %d dB, nsim = %d, Modulation Mode %d \n',SNR(snrIdx),nsim,nMod);
        end
    end
    
end
nmseSV = reshape(mean(NMSESV,2),[],modNum);
% figure(figCount);figCount = figCount + 1;
% semilogy(SNR,nmseSV);  
%% Baseline2
clear stateAL
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1e5;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 100;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
% inputParam.lambdaS = 0.005;
% inputParam.JESD = 1;
% inputParam.NiterA = 1;
% inputParam.lambda0 = 0.02;
% inputParam.convBreak = 0;
for nMod = 1:modNum
    inputParam.constell = constell(nMod).Cons;
    for snrIdx =1:length(SNR)
        inputParam.gammaOmega = 1e6*(mean(signalPower)/snr(snrIdx))^(-1);
        for nsim = 1:Nsim
            yp = zp + sqrt(signalPowerp/snr(snrIdx)).*noisedatap(:,:,nsim);
            yd = zd(:,:,nMod) + sqrt(signalPowerd/snr(snrIdx)).*noisedatad(:,:,nsim);
            [xestAL(:,:,snrIdx,nsim,nMod),destAL(:,:,snrIdx,nsim,nMod),stateAL(snrIdx,nMod)] = ALESD(yp,yd,Phip,Phid,inputParam);
            NMSEAL(snrIdx,nsim,nMod) = mean(vecnorm(xestAL(:,:,snrIdx,nsim,nMod) - x,2,1).^2./vecnorm(x,2,1).^2);
            dDecodedAL(:,:,snrIdx,nsim,nMod) = mleDecoder(destAL(:,:,snrIdx,nsim,nMod),constell(nMod).Cons);
            SERAL(snrIdx,nsim,nMod) = sum(dDecodedAL(:,:,snrIdx,nsim,nMod)~=d(:,:,nMod),'all')/numel(d(:,:,nMod));
            fprintf('Baseline 2, SNR = %d dB, nsim = %d, Modulation Mode %d.\n',SNR(snrIdx),nsim,nMod);
        end
    end
end
nmseAL = reshape(mean(NMSEAL,2),[],modNum);
% 
% SERAL  =  sum(destAL~=d,'all')/numel(d)

%% Baseline 3
clear stateWOD
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1e5;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 100;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
inputParam.JESD = 0;
% inputParam.NiterA = 1;
% inputParam.lambda0 = 0.02;
% inputParam.convBreak = 0;
% [xestWOD,destWOD,stateWOD] = MVVBI(yp,yd,Phip,Phid,inputParam);

for snrIdx =1:length(SNR)
    inputParam.gammaOmega = 1e6*(mean(signalPower)/snr(snrIdx))^(-1);
    for nsim = 1:Nsim
        yp = zp + sqrt(signalPowerp/snr(snrIdx)).*noisedatap(:,:,nsim);
        yd = zd(:,:,1) + sqrt(signalPowerd/snr(snrIdx)).*noisedatad(:,:,nsim);
        [xestWOD(:,:,snrIdx,nsim),destWOD(:,:,snrIdx,nsim),stateWOD(snrIdx)] = MVVBI(yp,yd,Phip,Phid,inputParam);
        
        NMSEWOD(snrIdx,nsim) = mean(vecnorm(xestWOD(:,:,snrIdx,nsim) - x,2,1).^2./vecnorm(x,2,1).^2);
        fprintf('Baseline 3, SNR = %d dB, nsim = %d \n',SNR(snrIdx),nsim);
    end
end
nmseWOD = reshape(mean(NMSEWOD,2),[],1);

%% Proposed Algorithm
dbstop if error;
clear stateProposed
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1e5;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 100;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
inputParam.JESD = 1;
% inputParam.NiterA = 1;
% inputParam.lambda0 = 0.02;
% inputParam.convBreak = 0;
for nMod = 1:modNum
    inputParam.constell = constell(nMod).Cons;
    for snrIdx =1:length(SNR)
        inputParam.gammaOmega = 1e6*(mean(signalPower)/snr(snrIdx))^(-1);
        for nsim = 1:Nsim
            yp = zp + sqrt(signalPowerp/snr(snrIdx)).*noisedatap(:,:,nsim);
            yd = zd(:,:,nMod) + sqrt(signalPowerd/snr(snrIdx)).*noisedatad(:,:,nsim);
            [xestProposed(:,:,snrIdx,nsim,nMod),destProposed(:,:,snrIdx,nsim,nMod),stateProposed(snrIdx,nMod)] = MVVBI(yp,yd,Phip,Phid,inputParam);
            NMSEProposed(snrIdx,nsim,nMod) = mean(vecnorm(xestProposed(:,:,snrIdx,nsim,nMod) - x,2,1).^2./vecnorm(x,2,1).^2);
            dDecodedProposed(:,:,snrIdx,nsim,nMod) = mleDecoder(destProposed(:,:,snrIdx,nsim,nMod),constell(nMod).Cons);
            SERProposed(snrIdx,nsim,nMod) = sum(dDecodedProposed(:,:,snrIdx,nsim,nMod)~=d(:,:,nMod),'all')/numel(d(:,:,nMod));
            fprintf('MVVBI, SNR = %d dB, nsim = %d, Modulation Mode %d. \n',SNR(snrIdx),nsim,nMod);
        end
    end
end
nmseProposed = reshape(mean(NMSEProposed,2),[],3);

% figure(figCount);figCount = figCount + 1;
% for tt = 1:Nt
%     stem(state.st(:,tt));hold on;
% end
% set(gca,'YScale','log');

% vecnorm(xest-gg,2).^2
% if inputParam.JESD
%     dDecoded = mleDecoder(dest,constell);
%     SERProposed  =  sum(dDecoded~=d,'all')/numel(d)
% end
% state.telap
%%
toc
modStr = strings(3,1);
modStr(1) = '4QAM';modStr(2) = '8PSK';modStr(3) = '16QAM';
%%
for nMod = 1:modNum
    figure(figCount);figCount = figCount + 1;
    semilogy(SNR,nmseWOD(:,1),'-.d');
    hold on;

    semilogy(SNR,nmseVBI(:,nMod),':s');
    semilogy(SNR,nmseSV(:,nMod),':>');
    semilogy(SNR,nmseAL(:,nMod),'--<');
    semilogy(SNR,nmseProposed(:,nMod),'-o');
    legend('MVSBLWOD','VBI','SVSBL','ALESD','MVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');
    ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
    xlabel('SNR',Interpreter='tex',FontName='Times New Roman',fontsize=15);
    grid on;
    nameString = strcat('fig\figureSNR',modStr(nMod),'.fig');
    savefig(nameString);
end

%%
for nMod = 1:modNum
    figure(figCount);figCount = figCount + 1;

    semilogy(SNR,mean(SERProposed(:,:,nMod),2),'-o');hold on;
    semilogy(SNR,mean(SERVBI(:,:,nMod),2),':s');
    
    semilogy(SNR,mean(SERSV(:,:,nMod),2),':s');
    semilogy(SNR,mean(SERAL(:,:,nMod),2),'--<');

    legend('MVSBL','VBI','SVSBL','ALESD',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');
    ylabel('SER',Interpreter='tex',FontName='Times New Roman',FontSize=15);
    xlabel('SNR',Interpreter='tex',FontName='Times New Roman',Fontsize=15);grid on;
    nameString = strcat('fig\figureSER',modStr(nMod),'.fig');
    savefig(nameString);
end
%%
datestr(now)
%% Auto save
str = datestr(now);
str = strrep(str,' ','');
str = strrep(str,'-','');
str = strrep(str,':','');
str = strcat('results\',str,'mainSNRtest','.mat');
save(str,'xestProposed','destProposed','destVBI','xestWOD','xestVBI','xestAL','destAL','xestSV',...
    'destSV','NMSEWOD','NMSEAL','NMSEVBI','NMSESV','NMSEProposed','SERAL','SERVBI','SERSV','SERProposed');