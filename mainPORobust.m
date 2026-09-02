%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Position offset\n');
% addpath('boundedline','singlepatch','catuneven','Inpaint_nans');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)
Ns = 1024;
deltaF = 30e3;   % Subcarrier spacing  
BW = Ns*deltaF;% BW = 50e6;
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
load('dataset\3\PhipGenie.mat');
load('dataset\3\PhiG');

% load('dataset\1\xUAVGT.mat');                                                % Ground Truth
load('dataset\1\xUAV.mat'); 
% load('dataset\1\tauGT.mat'); 


load("dataset\1\noisedatap.mat");
load("dataset\xOG.mat");
load("dataset\3\zp.mat");
load("dataset\1\derivP.mat");
load('dataset\3\derivG');

load("dataset\1\sMGParam.mat");
load('dataset\3\sMGParamG');

load('dataset\1\grid.mat');
load('dataset\1\xR.mat');
load('dataset\1\ULALine.mat');
load('dataset\1\gridPre.mat');

% sMGParam.derivD = derivD;
% sMGParam.derivP = derivP;
% sMGParamG.derivP = derivG;

% d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);
numBS = size(xR,2);
tauGT = 0.25/BW;
arrayPO = [2,4,6,8,10];
xing = [1,0,-1, 0;
        0,1, 0,-1;
        0,0, 0, 0];
for indPO = 1:length(arrayPO)
    for indXing = 1:4
        xUAVGT(:,indPO,indXing) = xUAV + xing(:,indXing)*arrayPO(indPO);
    end
end
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
inputParam.muPrior = zeros(size(x));
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
baselineSelect = [1,1,0,1,1]; %Proposed NoSL SV oracle  w/opp
% baselineSelect = [1,0,0,1,0]; %Proposed NoSL SV oracle  w/opp
% baselineSelect = [1,0,0,0]; %Proposed NoSL SV oracle

%%
tic
for nsim = 1:10
    inputParamtemp = inputParam;
    rmsePhaseI     = zeros(length(arrayPO),4);
    rmseProposed   = zeros(length(arrayPO),4);
    rmseNoSL       = zeros(length(arrayPO),4);
    rmseSV         = zeros(length(arrayPO),4);
    rmseWOPP       = zeros(length(arrayPO),4);
    rmseOracle     = zeros(length(arrayPO),4);
%     xProposed      = cell(length(tauGT),1);
%     xNoSL          = cell(length(tauGT),1);
%     xSV            = cell(length(tauGT),1);
%     xNoAid         = cell(length(tauGT),1);
%     xOracle        = cell(length(tauGT),1);
%     stateProposed  = cell(length(tauGT),1);
%     tauProposed    = cell(length(tauGT),1);
%     xUAVProposed   = cell(length(tauGT),1);
%     tauSV          = cell(length(tauGT),1);
%     xUAVSV         = cell(length(tauGT),1);
    PhipAP    = cell(length(arrayPO),4);

    for expIdx = 1 : length(arrayPO)
        for xingIdx = 1:4
    %       for expIdx = 1 : 1
            signalPowerp = mean(vecnorm(zpPO{expIdx,xingIdx},2).^2/size(zpPO{expIdx,xingIdx},1)); 
            yp = zpPO{expIdx,xingIdx} + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,nsim);
            inputParamtemp.signalPowerp = signalPowerp;
            inputParamtemp.gg = x;
            inputParamtemp.gammaPrior = 1e+3;
            inputParamtemp.gammaOmega = 1e+3;
            inputParamtemp.sMGParam = sMGParam;
            inputParamtemp.muPrior = 0;
            %% Phase I
            if any(baselineSelect([1,3]))
    %           if 1
                inputPhase1.snr = snr;
                inputPhase1.numAntennas = sMGParam.Na;
                inputPhase1.subcarrierIndices = sMGParam.subcarrInd;
                inputPhase1.subcarrierStride = median(diff(sMGParam.subcarrInd));
                inputPhase1.deltaF = deltaF;
                inputPhase1.oversampling = 8;
                inputPhase1.numTargets = nnz(x(:,1));
                inputPhase1.grid = gridPre{1};
                [centroids{1}] = SpecEstimation2D(yp(:,1),inputPhase1);
                theta{1} = pi/2 - acos(centroids{1}(:,1));
                range{1} = centroids{1}(:,2)/2/deltaF/inputPhase1.subcarrierStride*c0;
            
                inputPhase1.grid = gridPre{2};
                [centroids{2}] = SpecEstimation2D(yp(:,2),inputPhase1 );
                theta{2} = acos(centroids{2}(:,1));
                range{2} = centroids{2}(:,2)/2/deltaF/inputPhase1.subcarrierStride*c0;
                [xUAVEst, bias,final_pos] = estimate_bias_by_coherence(xUAV(1:2)', xR(1:2,1:2)', theta{1}, range{1}, theta{2}, range{2});
                rmsePhaseI(expIdx,xingIdx) = evaluate_offgrid_performance(final_pos.', ones(1,length(final_pos)), inputParamtemp.locGT.', 5);
                tauEst = bias/c0;
                [~, indGrid] = find_nearest_grid_manual([sMGParam.grid.x;sMGParam.grid.y], final_pos.');
                sMGParam1 = sMGParam;
    %             sMGParam1.grid.x(indGrid) = final_pos(:,1); 
    %             sMGParam1.grid.y(indGrid) = final_pos(:,2);
                indGrid = get_neighbor_indices(indGrid, Nsize);
                indGrid = sort(indGrid,'ascend');
                sMGParam1.grid.x = sMGParam1.grid.x(indGrid);
                sMGParam1.grid.y = sMGParam1.grid.y(indGrid);
                sMGParam1.xUAVinput(1:2,:) = xUAVEst';
                sMGParam1.vecTO = tauEst;
                
            
                for indBS = 1:numBS
                    sMGParam1.xR = xR(:,indBS);
                    sMGParam1.ULALine = ULALine(:,indBS);
                %                 sMGParam.symbol = ones(length(pilotSubcarrInd),1);
                    [PhipAP{expIdx,xingIdx}(:,:,indBS),sMGParam1.derivP(indBS)] = sensingMatrixGenWSymbolOG(sMGParam1);
                end
            else
                indGrid = 1:Nsize^2;
            end
    
    
           
            if baselineSelect(1)      
                inputParamtemp.gg = x(indGrid,:);
                inputParamtemp.sMGParam = sMGParam1;
                inputParamtemp.flagSL = 1;
                inputParamtemp.flagOG = 1;
                inputParamtemp.flagSV = 0;
        %     inputParam.gammaPrior = 1e+3;
    %             inputParamtemp.gammaOmega = 1e+3;
                [~,stateProposedtemp(expIdx,xingIdx)] = JSLES4(yp,PhipAP{expIdx,xingIdx},inputParamtemp);
                rmseProposed(expIdx,xingIdx) = stateProposedtemp(expIdx,xingIdx).RMSEfinal;
            end
    
            if baselineSelect(5)      
                inputParamtemp.gg = x(:,:);
                inputParamtemp.sMGParam = sMGParam;
                inputParamtemp.flagSL = 1;
                inputParamtemp.flagOG = 1;
                inputParamtemp.flagSV = 0;
        %     inputParam.gammaPrior = 1e+3;
    %             inputParamtemp.gammaOmega = 1e+3;
                [~,stateWOPP(expIdx,xingIdx)] = JSLES4(yp,Phip,inputParamtemp);
                rmseWOPP(expIdx,xingIdx) = stateWOPP(expIdx,xingIdx).RMSEfinal;
            end
            %% NoAid
    %         if baselineSelect(5)
    %             inputParamtemp.gg = x(:,:);
    %             inputParamtemp.sMGParam = sMGParam;
    %             inputParamtemp.flagSL = 1;
    %             inputParamtemp.flagOG = 1;
    %             inputParamtemp.flagSV = 1;
    %             [~,stateNoAidtemp(expIdx)] = JSLES4(yp(:,1),Phip(:,:,1),inputParamtemp);
    %             rmseNoAid(expIdx)    = stateNoAidtemp(expIdx).RMSEfinal;
    %         end
            %% OG-MVVAMP
    
            if baselineSelect(2)
                inputParamtemp.gg = x(:,:);
                inputParamtemp.sMGParam = sMGParam;
                inputParamtemp.flagSL = 0;
                inputParamtemp.flagOG = 1;
                inputParamtemp.flagSV = 0;
                [~,stateNoSLtemp(expIdx,xingIdx)] = JSLES4(yp,Phip,inputParamtemp);
                rmseNoSL(expIdx,xingIdx) = stateNoSLtemp(expIdx,xingIdx).RMSEfinal;
            end
    %         nmseNoSL{expIdx} = ...
    %             mean(vecnorm(xNoSL{expIdx} - x,2).^2./vecnorm(x,2).^2);
    %         NMSENoSL(:,nsim) = stateNoSL(nsim).NMSE;
        
            %% OG-SV
    
            if baselineSelect(3)
                inputParamtemp.gg = x(indGrid,:);
                inputParamtemp.sMGParam = sMGParam1;
                inputParamtemp.flagSL = 1;
                inputParamtemp.flagOG = 1;
                inputParamtemp.flagSV = 1;
    %             inputParamtemp.gammaOmega = 1e+2;
                [~,stateSVtemp(expIdx,xingIdx)] = JSLES4(yp,PhipAP{expIdx,xingIdx},inputParamtemp);
        %         tauSV{expIdx} = norm(statetmp.deltatau);
        %         xUAVSV{expIdx} = norm(statetmp.xUAV - xUAV,'fro');
                rmseSV(expIdx,xingIdx) = stateSVtemp(expIdx,xingIdx).RMSEfinal;
    %             inputParamtemp.gammaOmega = 1e+3;
    %             inputParamtemp.gammaOmega = (snr(expIdx) == 10)*1e3 + (snr(expIdx) < 10)*1e3;
            end
    %         nmseSV{expIdx} = ...
    %             mean(vecnorm(xSV{expIdx} - x,2).^2./vecnorm(x,2).^2);
    %         NMSESV(:,nsim) = stateSV(nsim).NMSE;
        
            %% Oracle
    
            if baselineSelect(4)
                inputParamtemp.gg = x(:,:);
                inputParamtemp.sMGParam = sMGParamGPO{expIdx,xingIdx};
                inputParamtemp.flagSL = 0;
                inputParamtemp.flagOG = 1;
                inputParamtemp.flagSV = 0;
    %             inputParamtemp.gammaOmega = 5e+2;
    
                [~,stateOracletemp(expIdx,xingIdx)] = JSLES4(yp,PhiGPO{expIdx,xingIdx},inputParamtemp);
                rmseOracle(expIdx,xingIdx) = stateOracletemp(expIdx,xingIdx).RMSEfinal;
            end


        end
    end
    if baselineSelect(1)
        RMSEProposed(:,:,nsim) = (rmseProposed);
%         XProposed{nsim} = xProposed; 
%         XNoAid{nsim} = xNoAid;
%         RMSENoAid(:,nsim)     = cell2mat(rmseNoAid);
    end
    if baselineSelect(5)
        RMSEWOPP(:,:,nsim) = (rmseWOPP);
%         XProposed{nsim} = xProposed; 
%         XNoAid{nsim} = xNoAid;
%         RMSENoAid(:,nsim)     = cell2mat(rmseNoAid);
    end
    if baselineSelect(2)
        RMSENoSL(:,:,nsim)     = (rmseNoSL);
%         XNoSL{nsim} = xNoSL;
    end
    if baselineSelect(3)
        RMSESV(:,:,nsim)       = (rmseSV);
%         XSV{nsim} = xSV;
    end
    if baselineSelect(4)
        RMSEOracle(:,:,nsim)   = (rmseOracle);
%         XOracle{nsim} = xOracle;
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
    semilogy(arrayPO,mean(mean(RMSENoSL(:,:),3),2),'-.>',LineWidth=1.5,Color = [138,43,226]/255);hold on;
    strLegend(baselineCounter) = 'OG-VAMP (w Fusion & Imperfect Tx prior)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(5)
    semilogy(arrayPO,mean(RMSEWOPP(:,:),2),'-+',LineWidth=1.5,Color = [250,128,114]/255);hold on;
    strLegend(baselineCounter) = 'JSSES (w/o Preprocessing)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(1)
%     semilogy(arrayPO,mean(RMSENoAid(:,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
    semilogy(arrayPO,mean(mean(RMSEProposed(:,:,[1:3,5:10]),3),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
%     strLegend(baselineCounter) = 'JSSES (No Aid)';
    strLegend(baselineCounter) = 'JSSES';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(4)
    semilogy(arrayPO,mean(mean(RMSEOracle(:,:),3),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
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
xlabel('Position offset (m)',Interpreter='latex',FontName='Times New Roman',FontSize=15);

nameStr = strcat('fig\figurePORMSE','.fig');
savefig(nameStr);


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

str = strcat('results\',str,'mainPORobust','.mat');
% save(str,'XProposed','XNoSL','XSV','XVAMP','XOracle',...
%     'NMSEProposed','NMSENoSL','NMSESV','NMSEVAMP','NMSEOracle',...
%     'RMSEtauProposed','RMSExUAVProposed','RMSExUAVSV','RMSEtauSV');
save(str,...
    'RMSEProposed','RMSEWOPP','RMSENoSL','RMSEOracle',"arrayPO");
fprintf('\n PO Simulation Finished!\n');
%
