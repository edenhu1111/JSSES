%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
clc;
clear;close all;
dbstop if error;
fprintf('Target number.\n');
% addpath('boundedline','singlepatch','catuneven','Inpaint_nans');
% global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)

Ns = 1024;
deltaF = 30e3;   % Subcarrier spacing  
BW = Ns*deltaF;
% T = 1 / deltaF;   % symbol duration
% Tcp = T / 4;      % cyclic prefix duration
% Ts = T + Tcp;     %   total symbol duration
% sigmaPos = 1/sqrt(3);
% sigma2 = 10^(-50/10);
xRange = 100;
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
load('dataset\1\PhipGenie.mat');
load('dataset\1\PhiG','PhiG');

load('dataset\1\xUAVGT.mat');                                                % Ground Truth
load('dataset\1\xUAV.mat'); 
load('dataset\1\tauGT.mat'); 


load("dataset\1\noisedatap.mat");
% load("dataset\xOG.mat");
% load("dataset\1\zp.mat");

load('dataset\1\sparArray','sparArray');
load('dataset\1\zpSpar','zp');
load('dataset\1\xSpar','x');

load("dataset\1\derivP.mat");
load('dataset\1\derivG','derivG');

load("dataset\1\sMGParam.mat");
load('dataset\1\sMGParamG','sMGParamG');

load('dataset\1\grid.mat');
load('dataset\1\xR.mat');
load('dataset\1\ULALine.mat');
load('dataset\1\gridPre.mat');


% sMGParam.derivD = derivD;
% sMGParam.derivP = derivP;
% sMGParamG.derivP = derivG;

% d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);
numBS = size(xR,2);
snr = 0; 
% nMod = 3;
 
% figCount = 1;
% signalPower = mean(vecnorm([zp;zd(:,:,nMod)],2).^2)/(size(zp,1) + size(zd,1));
% signalPowerd = mean(vecnorm(zd(:,:,nMod),2).^2)/size(zd,1);
inputParam.xR = xR;
inputParam.ULALine = ULALine;
inputParam.convBreaker = 1;
inputParam.xUAV = xUAV;

% inputParam.gammaPrior = 1e-1*ones(size(x));
inputParam.muPrior = 0;
inputParam.Lambda = 0.999;

inputParam.niter = 80;
inputParam.EMiter = 150;
inputParam.dampFacGam = 0.6;
inputParam.dampFac = 0.7;
inputParam.Normalization = 0;
inputParam.lambda0 = 0.90;
inputParam.lambdaS = 1e-3;

inputParam.armijoBeta1 = 0.5;
inputParam.armijoBeta2 = 0.5;
inputParam.armijoBeta3 = 0.5;
inputParam.iterArmijoMax = 10;
inputParam.armijoSigma = 0.0001;

inputParam.armijoRho   = 0.5;
inputParam.sMGParam = sMGParam;

inputParam.gammaPrior = 1e+3;
inputParam.gammaOmega = 1e+3;
% inputParam.signalPowerp = signalPowerp;
baselineSelect = [1,1,1,1,1]; %Proposed NoSL SV oracle
%%
tic
parfor nsim = 6:10
    inputParamtemp = inputParam;
    rmsePhaseI     = zeros(length(sparArray),1);
    rmseProposed   = zeros(length(sparArray),1);
    rmseNoSL       = zeros(length(sparArray),1);
    rmseSV         = zeros(length(sparArray),1);
    rmseNoAid      = zeros(length(sparArray),1);
    rmseOracle     = zeros(length(sparArray),1);
%     xProposed      = cell(length(snr),1);
%     xNoSL          = cell(length(snr),1);
%     xSV            = cell(length(snr),1);
%     xNoAid         = cell(length(snr),1);
%     xOracle        = cell(length(snr),1);

    stateProposedtemp      = cell(length(sparArray),1);
    stateNoSLtemp          = cell(length(sparArray),1);
    stateNoAidtemp         = cell(length(sparArray),1);
    stateSVtemp            = cell(length(sparArray),1);
    stateOracletemp        = cell(length(sparArray),1);
    PhipAP    = cell(length(sparArray),1);
    inputPhase1 = [];
    for expIdx = 1 : length(sparArray)
        signalPowerp = mean(vecnorm(zp{expIdx},2).^2/size(zp{expIdx},1)); 
        inputParamtemp.signalPowerp = signalPowerp;
        indTarget = find(x{expIdx}(:,1));
        inputParamtemp.locGT = [grid.x(indTarget);grid.y(indTarget)]';
        yp = zp{expIdx} + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,nsim);
        centroids = cell(numBS,1);
        theta     = cell(numBS,1);
        range     = cell(numBS,1);

%         inputParamtemp.gg = x;

%         inputParamtemp.sMGParam = sMGParam;
        %% Phase I
        if any(baselineSelect([1,3]))
%           if 1
            inputPhase1.snr = snr;
            inputPhase1.numAntennas = sMGParam.Na;
            inputPhase1.subcarrierIndices = sMGParam.subcarrInd;
            inputPhase1.subcarrierStride = median(diff(sMGParam.subcarrInd));
            inputPhase1.deltaF = deltaF;
            inputPhase1.oversampling = 8;
            inputPhase1.numTargets = nnz(x{expIdx}(:,1));
            inputPhase1.grid = gridPre{1};
            [centroids{1}] = SpecEstimation2D(yp(:,1),inputPhase1);
            theta{1} = pi/2 - acos(centroids{1}(:,1));
            range{1} = centroids{1}(:,2)/2/deltaF/inputPhase1.subcarrierStride*c0;
        
            inputPhase1.grid = gridPre{2};
            [centroids{2}] = SpecEstimation2D(yp(:,2),inputPhase1 );
            theta{2} = acos(centroids{2}(:,1));
            range{2} = centroids{2}(:,2)/2/deltaF/inputPhase1.subcarrierStride*c0;
            [xUAVEst, bias,final_pos] = estimate_bias_by_coherence(xUAV(1:2)', xR(1:2,1:2)', theta{1}, range{1}, theta{2}, range{2});
            rmsePhaseI(expIdx) = evaluate_offgrid_performance(final_pos.', ones(1,length(final_pos)), inputParamtemp.locGT.', 5);
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
                [PhipAP{expIdx}(:,:,indBS),sMGParam1.derivP(indBS)] = sensingMatrixGenWSymbolOG(sMGParam1);
            end
        else
            indGrid = 1:Nsize^2;
        end

        %% Proposed JVAMP

        if baselineSelect(1)      
            inputParamtemp.gg = x{expIdx}(indGrid,:);
            inputParamtemp.sMGParam = sMGParam1;
            inputParamtemp.flagSL = 1;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
    %     inputParam.gammaPrior = 1e+3;
%             inputParamtemp.gammaOmega = 1e+3;
            [~,stateProposedtemp{expIdx}] = JSLES4(yp,PhipAP{expIdx},inputParamtemp);
            rmseProposed(expIdx) = stateProposedtemp{expIdx}.RMSEfinal;
        end
        %% NoAid
        if baselineSelect(5)
            inputParamtemp.gg = x{expIdx}(:,:);
            inputParamtemp.sMGParam = sMGParam;
            inputParamtemp.flagSL = 1;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 1;
            [~,stateNoAidtemp{expIdx}] = JSLES4(yp(:,1),Phip(:,:,1),inputParamtemp);
            rmseNoAid(expIdx)    = stateNoAidtemp{expIdx}.RMSEfinal;
        end
        %% OG-VAMP

        if baselineSelect(2)
            inputParamtemp.gg = x{expIdx}(:,:);
            inputParamtemp.sMGParam = sMGParam;
            inputParamtemp.flagSL = 0;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
            [~,stateNoSLtemp{expIdx}] = JSLES4(yp,Phip,inputParamtemp);
            rmseNoSL(expIdx) = stateNoSLtemp{expIdx}.RMSEfinal;
        end
%         nmseNoSL{expIdx} = ...
%             mean(vecnorm(xNoSL{expIdx} - x,2).^2./vecnorm(x,2).^2);
%         NMSENoSL(:,nsim) = stateNoSL(nsim).NMSE;
    
        %% OG-SV

        if baselineSelect(3)
            inputParamtemp.gg = x{expIdx}(indGrid,:);
            inputParamtemp.sMGParam = sMGParam1;
            inputParamtemp.flagSL = 1;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 1;
%             inputParamtemp.gammaOmega = 1e+2;
            [~,stateSVtemp{expIdx}] = JSLES4(yp,PhipAP{expIdx},inputParamtemp);
    %         tauSV{expIdx} = norm(statetmp.deltatau);
    %         xUAVSV{expIdx} = norm(statetmp.xUAV - xUAV,'fro');
            rmseSV(expIdx) = stateSVtemp{expIdx}.RMSEfinal;
%             inputParamtemp.gammaOmega = 1e+3;
%             inputParamtemp.gammaOmega = (snr(expIdx) == 10)*1e3 + (snr(expIdx) < 10)*1e3;
        end
%         nmseSV{expIdx} = ...
%             mean(vecnorm(xSV{expIdx} - x,2).^2./vecnorm(x,2).^2);
%         NMSESV(:,nsim) = stateSV(nsim).NMSE;
    
        %% Oracle

        if baselineSelect(4)
            inputParamtemp.gg = x{expIdx}(:,:);
            inputParamtemp.sMGParam = sMGParamG;
            inputParamtemp.flagSL = 0;
            inputParamtemp.flagOG = 1;
            inputParamtemp.flagSV = 0;
%             inputParamtemp.gammaOmega = 5e+2;

            [~,stateOracletemp{expIdx}] = JSLES4(yp,PhiG,inputParamtemp);
            rmseOracle(expIdx) = stateOracletemp{expIdx}.RMSEfinal;
        end

    end
%%
    if any(baselineSelect([1,3]))
        RMSEPhaseI(:,nsim) = rmsePhaseI;
    end
    if baselineSelect(1)
        RMSEProposed(:,nsim) = rmseProposed;
%         XProposed{nsim} = xProposed; 
%         XNoAid{nsim} = xNoAid;
        RMSENoAid(:,nsim)     = rmseNoAid;
%         stateProposed{nsim} = stateProposedtemp;
%         stateNoAid{nsim} = stateNoAidtemp;
    end
    if baselineSelect(2)
        RMSENoSL(:,nsim)     = rmseNoSL;
%         XNoSL{nsim} = xNoSL;
%         stateNoSL{nsim} = stateNoSLtemp;
    end
    if baselineSelect(3)
        RMSESV(:,nsim)       = rmseSV;
%         XSV{nsim} = xSV;
%         stateSV{nsim} = stateSVtemp;
    end
    if baselineSelect(4)
        RMSEOracle(:,nsim)   = rmseOracle;
%         XOracle{nsim} = xOracle;
%         stateOracle{nsim} = stateOracletemp;
    end 
    %%
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
    semilogy(sparArray,mean(RMSENoSL(:,:),2),'-.>',LineWidth=1.5,Color = [138,43,226]/255);hold on;
    strLegend(baselineCounter) = 'OG-VAMP (w Fusion & Imperfect Tx prior)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(3)
    semilogy(sparArray,mean(RMSESV(:,:),2),'--<',LineWidth=1.5,Color = [255,128,0]/255);hold on;
    strLegend(baselineCounter) = 'JSSES (w/o Fusion)';
    baselineCounter = baselineCounter + 1;
end
if baselineSelect(1)
    semilogy(sparArray,mean(RMSENoAid(:,:),2),'-s',LineWidth=1.5,Color =  [220,20,60]/255);hold on;
    semilogy(sparArray,mean(RMSEProposed(:,:),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
    strLegend(baselineCounter) = 'JSSES (No Aid)';
    strLegend(baselineCounter+1) = 'JSSES';
    baselineCounter = baselineCounter + 2;
end
if baselineSelect(4)
    semilogy(sparArray,mean(RMSEOracle(:,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
    strLegend(baselineCounter) = 'Oracle';
    baselineCounter = baselineCounter + 1;
end

grid on;
ylim([3e-2,2]);
% strLegend(1) = 'VAMP (w Imperfect Tx prior)';
% ylim([2e-3,1]);
% semilogy(NMSE,'-o');
% hold on;

% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend(strLegend,Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='best');
ylabel('RMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Target number',Interpreter='tex',FontName='Times New Roman',FontSize=15);

nameStr = strcat('fig\figureTN','.fig');
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

str = strcat('results\',str,'mainTN','.mat');
% save(str,'XProposed','XNoSL','XSV','XVAMP','XOracle',...
%     'NMSEProposed','NMSENoSL','NMSESV','NMSEVAMP','NMSEOracle',...
%     'RMSEtauProposed','RMSExUAVProposed','RMSExUAVSV','RMSEtauSV');
save(str,...
    'RMSEProposed','RMSENoSL','RMSESV','RMSENoAid','RMSEOracle');
fprintf('\n Target Number Simulation Finished!\n');
%
