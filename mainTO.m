%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Time offset\n');
% addpath('boundedline','singlepatch','catuneven','Inpaint_nans');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)
BW = 1024*30e3;
% BW = 50e6;
% Ns = 1024;
% deltaF = BW/Ns;   % Subcarrier spacing  
% T = 1 / deltaF;   % symbol duration
% Tcp = T / 4;      % cyclic prefix duration
% Ts = T + Tcp;     %   total symbol duration
% sigmaPos = 1/sqrt(3);
% sigma2 = 10^(-50/10);
xRange = 80;
% niter = 20; 
Na = 32;
Nv = 1;
Nsize = 20; 
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
% xR = [-42; 0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = 3;
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
% pilotSubcarrInd = 0:32:Ns-1;
% dataSubcarrInd = [];

Nt = 1;
Nsim = 10;
indGrid = 1:Nsize^2;

% vecTO = 0/BW*(rand-1/2)*ones(Nr,1)/2; % offset is ignored
% vecTO = 1/BW*(rand-1/2);

% dArray = c0/fc/2;
% biasAngle = [0]';                                                         %degree
% biasAngle = [90]';                                                      %degree

% biasAngle = zeros(Nr,2); 
% vecTO = zeros(Nr,1);
% tauEST = zeros(size(vecTO));
% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
% xTarget = zeros(3,16);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation

load('dataset\1\Phip.mat');
load('dataset\2\PhipGenie.mat');
load('dataset\2\PhiG');

load('dataset\1\xUAVGT.mat');                                                % Ground Truth
load('dataset\1\xUAV.mat'); 
load('dataset\2\tauGT.mat'); 


load("dataset\1\noisedatap.mat");
load("dataset\xOG.mat");
load("dataset\2\zp.mat");
load("dataset\2\derivP.mat");
load('dataset\2\derivG');

load("dataset\1\sMGParam.mat");
load('dataset\2\sMGParamG');

load('dataset\1\grid.mat');
load('dataset\1\xR.mat');
load('dataset\1\ULALine.mat');

% sMGParam.derivD = derivD;
% sMGParam.derivP = derivP;
% sMGParamG.derivP = derivG;

% d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);
numBS = size(xR,2);
% tauGT = [-30,-20,-10,0,10,20,30]*1e-9;

snr = 0; 
% nMod = 3;
 
% figCount = 1;
% signalPower = mean(vecnorm([zp;zd(:,:,nMod)],2).^2)/(size(zp,1) + size(zd,1));
% signalPowerd = mean(vecnorm(zd(:,:,nMod),2).^2)/size(zd,1);
inputParam.xR = xR;
inputParam.ULALine = ULALine;
inputParam.convBreaker = 1;
inputParam.xUAV = xUAV;
indTarget = find(x(:,1));
inputParam.locGT = [grid.x(indTarget);grid.y(indTarget)]';
% inputParam.gammaPrior = 1e-1*ones(size(x));
inputParam.muPrior = 0;
inputParam.Lambda = 0.999;

inputParam.niter = 20;
inputParam.EMiter = 200;
inputParam.dampFacGam = 0.6;
inputParam.dampFac = 0.65;
inputParam.Normalization = 0;
inputParam.lambda0 = 0.90;
inputParam.lambdaS = 1e-4;

inputParam.armijoBeta1 = 0.5;
inputParam.armijoBeta2 = 0.5;
inputParam.armijoBeta3 = 0.5;
inputParam.iterArmijoMax = 10;
inputParam.armijoSigma = 0.0001;

inputParam.armijoRho   = 0.5;
inputParam.sMGParam = sMGParam;
baselineSelect = [1,1,0,1]; %Proposed NoSL SV oracle
% baselineSelect = [1,0,0,0]; %Proposed NoSL SV oracle

%%
tic
for nsim = 1:1
    inputParamtemp = inputParam;
    rmseProposed   = cell(length(tauGT),1);
    rmseNoSL       = cell(length(tauGT),1);
    rmseSV         = cell(length(tauGT),1);
    rmseNoAid      = cell(length(tauGT),1);
    rmseOracle     = cell(length(tauGT),1);
    xProposed      = cell(length(tauGT),1);
    xNoSL          = cell(length(tauGT),1);
    xSV            = cell(length(tauGT),1);
    xNoAid         = cell(length(tauGT),1);
    xOracle        = cell(length(tauGT),1);
    stateProposed  = cell(length(tauGT),1);
%     tauProposed    = cell(length(tauGT),1);
%     xUAVProposed   = cell(length(tauGT),1);
%     tauSV          = cell(length(tauGT),1);
%     xUAVSV         = cell(length(tauGT),1);


    for expIdx = 1 : length(tauGT)
%       for expIdx = 1 : 1
        signalPowerp = mean(vecnorm(zpTO{expIdx},2).^2/size(zpTO{expIdx},1)); 
        yp = zpTO{expIdx} + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,nsim);
        inputParamtemp.signalPowerp = signalPowerp;
        inputParamtemp.gg = x;
        inputParamtemp.gammaPrior = 1e+3;
        inputParamtemp.gammaOmega = 1e+3;
        inputParamtemp.sMGParam = sMGParam;
        inputParamtemp.muPrior = 0;

        %% Proposed JVAMP
%         inputParamtemp.sMGParam = sMGParam1;
        if baselineSelect(1)
            inputParamtemp.flagSL = 1;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
    %     
            [xProposed{expIdx},statetmp] = JSLES4(yp,Phip(:,indGrid,:),inputParamtemp);
            rmseProposed{expIdx} = statetmp.RMSEfinal;

%             inputParamtemp.flagSV = 1;
%             [xNoAid{expIdx},statetmp] = JSLES4(yp(:,1),Phip(:,indGrid,1),inputParamtemp);
%             rmseNoAid{expIdx} = statetmp.RMSEfinal;
%             inputParamtemp.gammaOmega = 1e+3;
        %% OG-MVVAMP
        end
        if baselineSelect(2)
            inputParamtemp.flagSL = 0;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
            [xNoSL{expIdx},statetmp] = JSLES4(yp,Phip(:,indGrid,:),inputParamtemp);
            rmseNoSL{expIdx} = statetmp.RMSEfinal;
        end
%         nmseNoSL{expIdx} = ...
%             mean(vecnorm(xNoSL{expIdx} - x,2).^2./vecnorm(x,2).^2);
%         NMSENoSL(:,nsim) = stateNoSL(nsim).NMSE;
    
        %% OG-SV
        if baselineSelect(3)
            inputParamtemp.flagSL = 1;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 1;
%             inputParamtemp.gammaOmega = 1e+2;
            [xSV{expIdx},statetmp] = JSLES4(yp,Phip(:,indGrid,:),inputParamtemp);
    %         tauSV{expIdx} = norm(statetmp.deltatau);
    %         xUAVSV{expIdx} = norm(statetmp.xUAV - xUAV,'fro');
            rmseSV{expIdx} = statetmp.RMSEfinal;
%             inputParamtemp.gammaOmega = 1e+3;
%             inputParamtemp.gammaOmega = (snr(expIdx) == 10)*1e3 + (snr(expIdx) < 10)*1e3;
        end
%         nmseSV{expIdx} = ...
%             mean(vecnorm(xSV{expIdx} - x,2).^2./vecnorm(x,2).^2);
%         NMSESV(:,nsim) = stateSV(nsim).NMSE;
    
        %% Oracle
        if baselineSelect(4)
            inputParamtemp.flagSL = 0;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
%             inputParamtemp.gammaOmega = 5e+2;
%             inputParamtemp.sMGParam = sMGParamG1;
            inputParamtemp.sMGParam = sMGParamGTO{expIdx};
%             inputParamtemp.sMGParam.derivP = derivGTO{expIdx};
            [xOracle{expIdx},statetmp] = JSLES4(yp,PhiGTO{expIdx}(:,indGrid,:),inputParamtemp);
            rmseOracle{expIdx} = statetmp.RMSEfinal;
        end

        %% VAMP
%         x1 = zeros(size(x));
%         for tt = 1:numBS
% %             inputParamtemp.gg = x(:,tt)*sqrt(10^((snr(expIdx)-10)/10));
%             inputParamtemp.gg = x(:,tt);
%             [x1(:,tt),~] = EMBGvampSVD(yp(:,tt),Phip(:,:,tt),inputParamtemp);
%         end
%         xVAMP{expIdx} = x1;
%         inputParamtemp.gg = x;
%         nmseVAMP{expIdx} = ...
%             mean(vecnorm(xVAMP{expIdx} - x,2).^2./vecnorm(x,2).^2);


    end
    if baselineSelect(1)
        RMSEProposed(:,nsim) = cell2mat(rmseProposed);
        XProposed{nsim} = xProposed; 
%         XNoAid{nsim} = xNoAid;
%         RMSENoAid(:,nsim)     = cell2mat(rmseNoAid);
    end
    if baselineSelect(2)
        RMSENoSL(:,nsim)     = cell2mat(rmseNoSL);
        XNoSL{nsim} = xNoSL;
    end
    if baselineSelect(3)
        RMSESV(:,nsim)       = cell2mat(rmseSV);
        XSV{nsim} = xSV;
    end
    if baselineSelect(4)
        RMSEOracle(:,nsim)   = cell2mat(rmseOracle);
        XOracle{nsim} = xOracle;
    end 
end
toc
%% Plot


%%
figure;
% semilogy(snr,mean(NMSEVAMP(:,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
% semilogy(snr,mean(NMSENoSL(:,:),2),'-.>',LineWidth=1.5,Color = [138,43,226]/255);hold on;
% semilogy(snr,mean(NMSESV(:,:),2),'--<',LineWidth=1.5,Color = [255,128,0]/255);hold on;
% semilogy(snr,mean(NMSENoAid(:,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
% semilogy(snr,mean(NMSEProposed(:,:),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
% semilogy(snr,mean(NMSEOracle(:,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);
strLegend = strings(5,1);baselineCounter = 1;
if baselineSelect(2)
    semilogy(tauGT'*1e9,mean(RMSENoSL(:,:),2),'-.>',LineWidth=1.5,Color = [138,43,226]/255);hold on;
    strLegend(baselineCounter) = 'OG-VAMP (w Fusion & Imperfect Tx prior)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(3)
    semilogy(tauGT'*1e9,mean(RMSESV(:,:),2),'--<',LineWidth=1.5,Color = [255,128,0]/255);hold on;
    strLegend(baselineCounter) = 'JSSES (w/o Fusion)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(1)
%     semilogy(tauGT'*1e9,mean(RMSENoAid(:,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
    semilogy(tauGT'*1e9,mean(RMSEProposed(:,:),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
%     strLegend(baselineCounter) = 'JSSES (No Aid)';
    strLegend(baselineCounter+1) = 'JSSES';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(4)
    semilogy(tauGT'*1e9,mean(RMSEOracle(:,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
    strLegend(baselineCounter) = 'Oracle';
    baselineCounter = baselineCounter + 1;
end

grid on;
ylim([5e-2,10]);
% strLegend(1) = 'VAMP (w Imperfect Tx prior)';
% ylim([2e-3,1]);
% semilogy(NMSE,'-o');
% hold on;

% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend(strLegend,Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='best');
ylabel('RMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Time offset',Interpreter='latex',FontName='Times New Roman',FontSize=15);

% nameStr = strcat('fig\figureTORMSE','.fig');
% savefig(nameStr);


%% borderline (abandoned)
% figure;
% set(gca, 'YScale', 'log');  % 将当前坐标轴的y轴设置为对数刻度
% 
% if baselineSelect(2)
% % bNoSL = [sqrt(mean((RMSENoSL - mean(RMSENoSL,2)).^2,2))];
% bNoSL = [-min((RMSENoSL - mean(RMSENoSL,2)),[],2),max((RMSENoSL - mean(RMSENoSL,2)),[],2)];
% boundedline(snr,mean(RMSENoSL,2)',bNoSL,'cmap',[138,43,226]/255);hold on;
% end
% if baselineSelect(3)
% bSV = [-min((RMSESV - mean(RMSESV,2)),[],2),max((RMSESV - mean(RMSESV,2)),[],2)];
% boundedline(snr,mean(RMSESV,2)',bSV,'cmap',[255,128,0]/255);hold on;
% end
% if baselineSelect(1)
% bProposed = [-min((RMSEProposed - mean(RMSEProposed,2)),[],2),max((RMSEProposed - mean(RMSEProposed,2)),[],2)];
% bNoAid = [-min((RMSENoAid - mean(RMSENoAid,2)),[],2),max((RMSENoAid - mean(RMSENoAid,2)),[],2)];
% boundedline(snr,mean(RMSENoAid,2)',bNoAid,'cmap',[220,20,60]/255);hold on;
% boundedline(snr,mean(RMSEProposed,2)',bProposed,'cmap',[0,0,255]/255);hold on;
% end
% if baselineSelect(4)
% bOracle = [-min((RMSEOracle - mean(RMSEOracle,2)),[],2),max((RMSEOracle - mean(RMSEOracle,2)),[],2)];
% boundedline(snr,mean(RMSEOracle,2)',bOracle,'cmap',[0,0,0]/255);hold on;
% % semilogy(snr,mean(RMSEOracle(:,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
% end
% 
% nameStr = strcat('fig\figureSNRRMSE','.fig');
% savefig(nameStr);
%%
str = datestr(now);
str = strrep(str,' ','');
str = strrep(str,'-','');
str = strrep(str,':','');

str = strcat('results\',str,'mainTO','.mat');
% save(str,'XProposed','XNoSL','XSV','XVAMP','XOracle',...
%     'NMSEProposed','NMSENoSL','NMSESV','NMSEVAMP','NMSEOracle',...
%     'RMSEtauProposed','RMSExUAVProposed','RMSExUAVSV','RMSEtauSV');
save(str,'XProposed','XNoSL','XOracle',...
    'RMSEProposed','RMSENoSL','RMSEOracle',"tauGT");
fprintf('\n TO Simulation Finished!\n');
%