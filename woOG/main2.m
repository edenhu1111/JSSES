%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Convergence Behavior\n');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle gg xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)

deltaF = 240e+3;   % subcarrier spacing
Ns = 450;
BW = deltaF*Ns;  
T = 1 / deltaF;   % symbol duration
Tcp = T / 4;      % cyclic prefix duration
Ts = T + Tcp;     % total symbol duration
sigmaPos = 1/sqrt(3);
sigma2 = 10^(-50/10);
xRange = 50;
niter = 20;
Na = 16;
Nv = 1;
Nsize = 20;
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
xR = [-22; -0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = size(xR,2);
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
ind = 0:Ns-1;
pilotSubcarrInd = 0:15:Ns-1;
dataSubcarrInd = 0:1:Ns-1;
dataSubcarrInd = setdiff(dataSubcarrInd,pilotSubcarrInd);
% dataSubcarrInd = pilotSubcarrInd;
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
Nt = 10;
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

load('dataset\1\Phid.mat');
load('dataset\1\Phip.mat');
load('dataset\1\xUAV.mat');

load("dataset\1\d.mat");
load("dataset\1\noisedatap.mat");
load("dataset\1\noisedatad.mat");
load("dataset\x.mat");
load("dataset\1\zp.mat");
load("dataset\1\zd.mat");

d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);

gg = x;

% figCount = 1;
signalPower = mean(vecnorm([zp;zd(:,:,1)]).^2)/(size(zp,1) + size(zd,1));
signalPowerp = mean(vecnorm(zp).^2)/size(zp,1);
signalPowerd = mean(vecnorm(zd(:,:,1)).^2)/size(zd,1);

snr = 10; nMod = 3;
yp = zp + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,1);
yd = zd(:,:,nMod) + sqrt(signalPowerd/(10^(snr/10)))*noisedatad(:,:,1);
%% VBI
clear stateSBL
figCount = 1;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
inputParam.constell = constell(nMod).Cons;
inputParam.gammaOmega = 1e6;
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
inputParam.convBreak = 0;

[xestVBI,destVBI,stateVBI] = VBI2(yp,yd,Phip,Phid,inputParam);
nmseVBI = reshape(stateVBI.NMSE,[],1);
nmseVBI = nmseVBI(nmseVBI~=0);
% NMSEVBI = vec(stateVBI.NMSE(end,:));
NMSEVBI = find4Conv(stateVBI.NMSE);
% NMSEVBI = [nmseVBI(1),NMSEVBI];
dDecodedVBI = mleDecoder(destVBI,constell(nMod).Cons);
SERVBI = sum(dDecodedVBI~=d(:,:,nMod),'all')/numel(d(:,:,nMod))


%% Baseline1 (SVSBL)
% figCount = 1;
% inputParam.constell = constell;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1e5;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
inputParam.Niter = 10;
inputParam.NiterA = 20;
inputParam.JESD = 1;
inputParam.SV = 1;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
% inputParam.NiterA = 20;
% inputParam.lambda0 = 0.5;
% inputParam.lambdaS = 0.002;

% inputParam.lambda0 = 0.8;
% inputParam.convBreak = 0;
% for tt = 1:Nt
%     gg = x(:,tt);
%     [xestSV(:,tt),destSV(:,tt),stateSV(tt)] = SVVBI(yp(:,tt),yd(:,tt),Phip(:,:,tt),Phid(:,:,tt),inputParam);
% end
[xestSV,destSV,stateSV] = SVVBI(yp,yd,Phip,Phid,inputParam);

nmseSV = reshape(stateSV.NMSE,[],1);
nmseSV = nmseSV(nmseSV~=0);
% NMSESV = vec(stateSV.NMSE(end,:));
NMSESV = find4Conv(stateSV.NMSE);
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
    SERSV  =  sum(dDecodedSV~=d(:,:,nMod),'all')/numel(d(:,:,nMod))
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
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1./sigma2;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 10;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
inputParam.JESD = 1;
% inputParam.NiterA = 10;
% inputParam.lambda0 = 0.05;
% inputParam.convBreak = 1;

[xestAL,destAL,stateAL] = ALESD(yp,yd,Phip,Phid,inputParam);
% figure(figCount);figCount = figCount + 1;
% for ii =1:Nt
%     subplot(1,Nt,ii);
%     stem(abs(xestAL(:,ii)),'o');
%     hold on;stem(abs(x(:,ii)),'x');
% end
% figure(figCount);figCount = figCount + 1;
nmseAL = reshape(stateAL.NMSE,[],1);
nmseAL = nmseAL(nmseAL~=0);
% NMSEAL = vec(stateAL.NMSE(end,:));
NMSEAL = find4Conv(stateAL.NMSE);
% semilogy(nmseAL);

SERAL  =  sum(destAL~=d(:,:,nMod),'all')/numel(d(:,:,nMod))

%% Baseline 3
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1./sigma2;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 10;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
inputParam.JESD = 0;
% inputParam.NiterA = 10;
% inputParam.lambda0 = 0.05;
% inputParam.convBreak = 1;
[xestWOD,destWOD,stateWOD] = MVVBI(yp,yd,Phip,Phid,inputParam);
% figure(figCount);figCount = figCount + 1;
% for ii =1:Nt
%     subplot(1,Nt,ii);
%     stem(abs(xestWOD(:,ii)),'o');
%     hold on;stem(abs(x(:,ii)),'x');
% end
% figure(figCount);figCount = figCount + 1;
nmseWOD = reshape(stateWOD.NMSE,[],1);
nmseWOD = nmseWOD(nmseWOD~=0);
% NMSEWOD = vec(stateWOD.NMSE(end,:));
NMSEWOD = find4Conv(stateWOD.NMSE);
% semilogy(nmseWOD);

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



%% Proposed Algorithm
dbstop if error;
gg = x;
inputParam.pilotSubcarrInd = pilotSubcarrInd;
inputParam.dataSubcarrInd  = dataSubcarrInd;
% inputParam.gammaOmega = 1./sigma2;
inputParam.dataNum = dataNum;
inputParam.eyeS = eyeS;
% inputParam.Niter = 10;
inputParam.xUAV = xUAV;
inputParam.xR = xR;
inputParam.JESD = 1;
% inputParam.NiterA = 10;
% inputParam.lambda0 = 0.05;
% inputParam.convBreak = 1;
[xest,dest,state] = MVVBI(yp,yd,Phip,Phid,inputParam);
% figure(figCount);figCount = figCount + 1;
% for ii =1:Nt
%     subplot(1,Nt,ii);
%     stem(abs(xest(:,ii)),'o');
%     hold on;stem(abs(x(:,ii)),'x');
% end
% figure(figCount);figCount = figCount + 1;
% nmse = reshape(state.NMSE,[],1);
% nmse = nmse(nmse~=0);
nmse = reshape(state.NMSE,[],1);
nmse = nmse(nmse~=0);
% NMSE = vec(state.NMSE(end,:));
NMSE = find4Conv(state.NMSE);
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
    dDecoded = mleDecoder(dest,constell(nMod).Cons);
    SERProposed  =  sum(dDecoded~=d(:,:,nMod),'all')/numel(d(:,:,nMod))
end
%%
modStr = strings(3,1);
modStr(1) = '4QAM';modStr(2) = '8PSK';modStr(3) = '16QAM';

figure(figCount);figCount = figCount + 1;

semilogy(nmse,'-o');hold on;semilogy(nmseVBI,':>');semilogy(nmseSV,':s');semilogy(nmseAL,'--<');semilogy(nmseWOD,'-.d');xlim([1 100]);

legend('MVSBL','VBI','SVVBL','ALESD','MVSBLWOD',Interpreter='tex',FontName='Times New Roman',FontSize=10);
ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
grid on;
nameStr = strcat('fig\figureConvergenceFull',modStr(nMod),'.fig');
savefig(nameStr);


figure(figCount);figCount = figCount + 1;

semilogy(NMSE,'-o');hold on;
semilogy(NMSEVBI,':>');
semilogy(NMSESV,':s');semilogy(NMSEAL,'--<');semilogy(NMSEWOD,'-.d');

legend('MVSBL','VBI','SVVBL','ALESD','MVSBLWOD',Interpreter='tex',FontName='Times New Roman',FontSize=10);xlim([1 10]);
ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
grid on;
nameStr = strcat('fig\figureConvergence',modStr(nMod),'.fig');
savefig(nameStr);
%%
str = datestr(now);
str = strrep(str,' ','');
str = strrep(str,'-','');
str = strrep(str,':','');

str = strcat('results\',str,'mainConvtest',modStr(nMod),'.mat');
save(str,'xest','dest','xestWOD','xestAL','destAL','xestSV',...
    'destSV','NMSEWOD','NMSEAL','NMSESV','NMSE','nmseWOD','nmseAL','nmseSV','nmse');

%