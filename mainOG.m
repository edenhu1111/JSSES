%% tests of MVSBI and baselines
% author: Yunbo Hu (Eden Hu)
%% Initalize the environment
% clc;
clear;close all;
dbstop if error;
fprintf('Convergence Behavior\n');
global c0 fc Na Nv deltaF sigma2 Nsize dArray biasAngle xRange
%% ISAC Tx & Rx
% System parameters
c0 = physconst('lightspeed');  % velocity of light
fc = 28e+9;                    % carrier frequency (Ka)

deltaF = 30e3;
% BW = 100e6;
Ns = 1024;
BW = Ns*deltaF;
% deltaF = BW/Ns;   % Subcarrier spacing  
% T = 1 / deltaF;   % symbol duration
% Tcp = T / 4;      % cyclic prefix duration
% Ts = T + Tcp;     %   total symbol duration
% sigmaPos = 1/sqrt(3);
% sigma2 = 10^(-50/10);
% xRange = 100;
niter = 20; 
% Na = 16;
% Nv = 1;
% xR = [-50,50,-0,0 ; -0,0,50,-50 ; 5,5,5,5];
% xR = [-42; 0; 5];

% xR = [-20,; -0; 5];
% xR = [-80,80,-0,0 ; -0,-0,80,-80 ; 2,2,2,2];

Nr = 1;
% subcarrInd = [0:2:63,64:6:191,192:4:255];
% subcarrInd = sort(randperm(256,64));
% pilotSubcarrInd = 0:32:Ns-1;
% dataSubcarrInd = [];



Nt = 2;
Nsim = 10;

% vecTO = 0/BW*(rand-1/2)*ones(Nr,1)/2; % offset is ignored
% vecTO = 1/BW*(rand-1/2);
                                                     %degree
% biasAngle = [90]';                                                      %degree

% biasAngle = zeros(Nr,2); 
% vecTO = zeros(Nr,1);
% tauEST = zeros(size(bias));
% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
xTarget = zeros(3,16);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation

load('dataset\1\Phip.mat');
load('dataset\1\PhipGenie.mat');
load('dataset\1\PhiG','PhiG');

load('dataset\1\xUAVGT.mat');                                                % Ground Truth
load('dataset\1\xUAV.mat'); 
load('dataset\1\tauGT.mat'); 


load("dataset\1\noisedatap.mat");
load("dataset\xOG.mat");
load("dataset\1\zp.mat");
% load("dataset\1\derivP.mat");
% load('dataset\1\derivG','derivG');

load("dataset\1\sMGParam.mat");
load('dataset\1\sMGParamG','sMGParamG');

load('dataset\1\grid.mat');
load('dataset\1\xR.mat');
load('dataset\1\ULALine.mat');
load('dataset/1/gridPre.mat');
xRange = 80;
Nsize = 20;
% d(:,:,2) = mleDecoder(d(:,:,2),constell(2).Cons);
numBS = size(xR,2);
snr = 0;  % UNIT: dB
% nMod = 3;
% dd = d(:,:,nMod);

% figCount = 1;
% signalPower = mean(vecnorm([zp;zd(:,:,nMod)],2).^2)/(size(zp,1) + size(zd,1));
signalPowerp = mean(vecnorm(zp,2).^2)/size(zp,1); 
% signalPowerd = mean(vecnorm(zd(:,:,nMod),2).^2)/size(zd,1);
inputParam.xR = xR;
inputParam.ULALine = ULALine;
inputParam.convBreaker = 0;
inputParam.xUAV = xUAV;
indTarget = find(x(:,1));
inputParam.locGT = [grid.x(indTarget);grid.y(indTarget)]';

% baselineVec = [1,1,1,1,0];% oracle, jesd, OG-VAMP, jesd(no fusion),VAMP
baselineVec = [1,1,1,1,0];% oracle, jesd, OG-VAMP, jesd(no fusion),VAMP

% snr = 10; nMod = 3;
tic
indGrid = 1:Nsize^2;
fprintf('SNR = %d dB.\n',snr);
xProposed = zeros(Nsize^2,numBS,Nsim);
xOracle = zeros(Nsize^2,numBS,Nsim);
xNoAid = zeros(Nsize^2,1,Nsim);
xSV = zeros(Nsize^2,numBS,Nsim);
xNoSL = zeros(Nsize^2,numBS,Nsim);
x1 = zeros(Nsize^2,numBS,Nsim);
inputParam.gg = x;

%% 
% dRange(1,:) = calc_signal_dist_range(xUAV(1:2).',xR(1:2,1)',[grid.x.',grid.y.']);
% dRange(2,:) = calc_signal_dist_range(xUAV(1:2).',xR(1:2,2)',[grid.x.',grid.y.']);
% ddRange(1,:) = dRange(1,1):c0/BW/4:dRange(1,2);
% ddRange(2,:) = dRange(1,1):c0/BW/4:dRange(1,2);
% 
% gridtau(1,:) = 2*deltaF*32*ddRange(1,:)/c0;
% gridtau(2,:) = 2*deltaF*32*ddRange(2,:)/c0;
% 
% gridangle(1,:) = linspace(-sqrt(3)/2,sqrt(3)/2,64);
% gridangle(2,:) = linspace(-sqrt(3)/2,sqrt(3)/2,64);
% 
% 
% gridPre{1}(1,:) = kron(gridangle(1,:),ones(size(gridtau(1,:))));
% gridPre{1}(2,:) = kron(ones(size(gridangle(1,:))),gridtau(1,:));
% 
% gridPre{2}(1,:) = kron(gridangle(2,:),ones(size(gridtau(2,:))));
% gridPre{2}(2,:) = kron(ones(size(gridangle(2,:))),gridtau(2,:));
%%
for nsim = 1:5
    yp = zp + sqrt(signalPowerp/(10^(snr/10)))*noisedatap(:,:,nsim);
    %% Phase1: tau_estimation
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
    [centroids{2}] = SpecEstimation2D(yp(:,2),inputPhase1);
    theta{2} = acos(centroids{2}(:,1));
    range{2} = centroids{2}(:,2)/2/deltaF/inputPhase1.subcarrierStride*c0;
    [xUAVEst, bias,final_pos] = estimate_bias_by_coherence(xUAV(1:2)', xR(1:2,1:2)', theta{1}, range{1}, theta{2}, range{2});
    tauEst = bias/c0;
    [~, indGrid] = find_nearest_grid_manual([sMGParam.grid.x;sMGParam.grid.y], final_pos.');
    sMGParam1 = sMGParam;
%     sMGParam1.grid.x(indGrid) = final_pos(:,1); 
%     sMGParam1.grid.y(indGrid) = final_pos(:,2);
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
        [PhipAP{nsim}(:,:,indBS),sMGParam1.derivP(indBS)] = sensingMatrixGenWSymbolOG(sMGParam1);
        

    end
%     yd = zd(:,:,nMod) + sqrt(signalPowerd/(10^(snr/10)))*noisedatad(:,:,nsim);
    %% Initialization
    inputParam.gammaPrior = 1e3;
    inputParam.muPrior = 0;
    inputParam.Lambda = 0.999;
    inputParam.gammaOmega = 1e+3;
    inputParam.niter = 50;
    Niter = 100;
    inputParam.EMiter = Niter;
    inputParam.dampFacGam = 0.6;
    inputParam.dampFac = 0.7;
    inputParam.Normalization = 0;
    inputParam.lambda0 = 0.90;
    inputParam.lambdaS = 0.001;

    inputParam.armijoBeta1 = 0.5;
    inputParam.armijoBeta2 = 0.5;
    inputParam.armijoBeta3 = 0.5;
    inputParam.iterArmijoMax = 10;
    inputParam.armijoSigma = 0;
 
    inputParam.armijoRho   = 0.5;
    inputParam.sMGParam = sMGParamG;
    %%
    %% Grid Initialization
%     [outTheta1] = CoarseEstimation(yp(:,1),inputParam);
%     outTheta1 = -outTheta1 + pi/2;
%     [outTheta2] = CoarseEstimation(yp(:,2),inputParam);
%     outTheta2 = outTheta2;
%     targets = calculate_possible_targets(xR(1:2,1),0,outTheta1 , xR(1:2,2), 0, outTheta2);
%     valid_idx = all(abs(targets) <= xRange/2,1);
%     targets = targets(:,valid_idx);
%     [~, indGrid] = find_nearest_grid_manual([sMGParam.grid.x;sMGParam.grid.y], targets);
% %     indGrid = indTarget;
%     sMGParamG1 = sMGParamG; sMGParam1 = sMGParam; 
%     sMGParamG1.grid.x = sMGParamG1.grid.x(indGrid); sMGParamG1.grid.y = sMGParamG1.grid.y(indGrid);
%     sMGParam1.grid.x = sMGParam1.grid.x(indGrid); sMGParam1.grid.y = sMGParam1.grid.y(indGrid);
%     for ii = 1:numBS
% 
%         sMGParamG1.derivP(ii).Phidx = sMGParamG1.derivP(ii).Phidx(:,indGrid,:);
%         sMGParamG1.derivP(ii).Phidy = sMGParamG1.derivP(ii).Phidy(:,indGrid,:);
%         sMGParamG1.derivP(ii).PhidxUAV = sMGParamG1.derivP(ii).PhidxUAV(:,indGrid,:);
%         sMGParamG1.derivP(ii).PhidyUAV = sMGParamG1.derivP(ii).PhidyUAV(:,indGrid,:);
%         sMGParamG1.derivP(ii).PhidzUAV = sMGParamG1.derivP(ii).PhidzUAV(:,indGrid,:);
%         sMGParamG1.derivP(ii).Phidtau  = sMGParamG1.derivP(ii).Phidtau(:,indGrid,:);
% 
% 
%         sMGParam1.derivP(ii).Phidx = sMGParam1.derivP(ii).Phidx(:,indGrid,:);
%         sMGParam1.derivP(ii).Phidy = sMGParam1.derivP(ii).Phidy(:,indGrid,:);
%         sMGParam1.derivP(ii).PhidxUAV = sMGParam1.derivP(ii).PhidxUAV(:,indGrid,:);
%         sMGParam1.derivP(ii).PhidyUAV = sMGParam1.derivP(ii).PhidyUAV(:,indGrid,:);
%         sMGParam1.derivP(ii).PhidzUAV = sMGParam1.derivP(ii).PhidzUAV(:,indGrid,:);
%         sMGParam1.derivP(ii).Phidtau  = sMGParam1.derivP(ii).Phidtau(:,indGrid,:);
%     end
%     inputParam.gammaPrior = inputParam.gammaPrior(indGrid,:);
%     inputParam.muPrior = inputParam.muPrior(indGrid,:);
%     inputParam.gg = x(indGrid,:);


    %% Oracle
%     inputParam.sMGParam = sMGParamG1;
    inputParam.sMGParam = sMGParamG;
    inputParam.gg = x;
%     indGrid = 1:Nsize^2;
    if baselineVec(1)
        inputParam.flagSL = 0;
        inputParam.flagOG = 1;
        inputParam.flagSV = 0;
        [xOracle(:,:,nsim),stateOracle(nsim)] = JSLES4(yp,PhiG(:,:,:),inputParam);
%         [xOracle(indGrid,:,nsim),stateOracle(nsim)] = JSLES2(yp(:,1),PhiG(:,indGrid,1),inputParam);
%         [xOracle(:,:,nsim),stateOracle(nsim)] = JSLES2(yp,PhipGenie,inputParam);
        NMSEOracle(:,nsim) = stateOracle(nsim).NMSEz;
        RMSEOracle(:,nsim) = stateOracle(nsim).RMSE;
%         figure;
%         scatter(stateOracle(nsim).xGrid,stateOracle(nsim).yGrid);hold on;scatter(grid.x,grid.y);
    end
    %%
%     figure;
%     surf(abs(vec2imag(xOracle(:,1,10))));shading interp;
    %% Proposed
%     inputParam.sMGParam = sMGParam1;
    inputParam.sMGParam = sMGParam1;
    inputParam.gg = x(indGrid,:);
    if baselineVec(2)
        inputParam.flagSL = 1;
        inputParam.flagOG = 1;
        inputParam.flagSV = 0; 
        [xProposed(indGrid,:,nsim),stateProposed(nsim)] = JSLES4(yp,PhipAP{nsim}(:,:,:),inputParam);
        NMSEProposed(:,nsim) = stateProposed(nsim).NMSEz;
        RMSEProposed(:,nsim) = stateProposed(nsim).RMSE;

        inputParam.gg = x;
        inputParam.sMGParam = sMGParam; 
        inputParam.flagSL = 1;
        inputParam.flagOG = 1;
        inputParam.flagSV = 1;
        [xNoAid(:,:,nsim),stateNoAid(nsim)] = JSLES4(yp(:,1),Phip(:,:,1),inputParam);
        NMSENoAid(:,nsim) = stateNoAid(nsim).NMSEz;
        RMSENoAid(:,nsim) = stateNoAid(nsim).RMSE;
    end
    %% OG-MVVAMP
    inputParam.gg = x;
    inputParam.sMGParam = sMGParam;
    if baselineVec(3)
        inputParam.flagSL = 0;
        inputParam.flagOG = 1;
        inputParam.flagSV = 0;
    
        [xNoSL(:,:,nsim),stateNoSL(nsim)] = JSLES4(yp,Phip(:,:,:),inputParam);
        NMSENoSL(:,nsim) = stateNoSL(nsim).NMSEz;
        RMSENoSL(:,nsim) = stateNoSL(nsim).RMSE;
    end

    %% OG-SVVAMP
    inputParam.gg = x(indGrid,:);
    inputParam.sMGParam = sMGParam1;
    if baselineVec(4)
        inputParam.flagSL = 1;
        inputParam.flagOG = 1;
        inputParam.flagSV = 1;
        [xSV(indGrid,:,nsim),stateSV(nsim)] = JSLES4(yp,PhipAP{nsim}(:,:,:),inputParam);
        NMSESV(:,nsim) = stateSV(nsim).NMSEz;
        RMSESV(:,nsim) = stateSV(nsim).RMSE;
    end

    %% VAMP
    if baselineVec(5)
        inputParam.flagSL = 0;
        inputParam.flagOG = 0;
        inputParam.flagSV = 1;
        for tt = 1:numBS
            inputParam.gg = x(:,tt);
            [x1(:,tt,nsim),statetmp] = EMBGvampSVD(yp(:,tt),PhipAP{nsim}(:,:,tt),inputParam);
    %         [x1(:,tt,nsim),statetmp] = EMBGvampSVD(yp(:,tt),PhipGenie(:,:,tt),inputParam);
            nmseVAMP(:,nsim,tt) = statetmp.NMSE;
        end
        NMSEVAMP(:,nsim) = mean(nmseVAMP(:,nsim,:),3);
        inputParam.gg = x;
    end

end
toc
%% Plot
% figure;
% surf(abs(vec2imag(xProposed(:,1))));shading interp;
% figure;stem(abs(xProposed(:,1)));hold on; stem(abs(x(:,1)));
if snr == 0
    gridxInit = sMGParam.grid.x;
    gridyInit = sMGParam.grid.y;

    gridxInit(indGrid) = stateProposed(1).xGrid;
    gridyInit(indGrid) = stateProposed(1).yGrid;
%     figure;scatter(stateProposed(1).xGrid,stateProposed(1).yGrid,'ks'); hold on;scatter(grid.x(indTarget),grid.y(indTarget),'r+');
    figure;scatter(gridxInit,gridyInit,'ks'); hold on;scatter(grid.x(indTarget),grid.y(indTarget),'r+');

    % % % scatter(sMGParam.grid.x(indTarget),sMGParam.grid.y(indTarget),'bs');
    legend('Grid location','Ground truth',Interpreter='latex',FontSize=15);xlim([-40,40]),ylim([-40,40]);xlabel('$x$ (m)',Interpreter='latex',FontSize=16);ylabel('$y$ (m)',Interpreter='latex',FontSize=16);
    savefig('fig/figureLoc.fig');
    print('-depsc2', 'fig/figureLoc.eps');
    
%     figure;scatter(stateSV(1).xGrid,stateSV(1).yGrid,'ks'); hold on;scatter(grid.x(indTarget),grid.y(indTarget),'r+');
%     % % % scatter(sMGParam.grid.x(indTarget),sMGParam.grid.y(indTarget),'bs');
%     legend('Grid location','Ground truth',Interpreter='latex',FontSize=15);xlim([-40,40]),ylim([-40,40]);xlabel('$x$ (m)',Interpreter='latex',FontSize=16);ylabel('$y$ (m)',Interpreter='latex',FontSize=16);
%     savefig('fig/figureLocSV.fig');
%     print('-depsc2', 'fig/figureLocSV.eps');
    
    figure;scatter(stateNoAid(1).xGrid,stateNoAid(1).yGrid,'ks'); hold on;scatter(grid.x(indTarget),grid.y(indTarget),'r+');
    % % % scatter(sMGParam.grid.x(indTarget),sMGParam.grid.y(indTarget),'bs');
    legend('Grid location','Ground truth',Interpreter='latex',FontSize=15);xlim([-40,40]),ylim([-40,40]);xlabel('$x$ (m)',Interpreter='latex',FontSize=16);ylabel('$y$ (m)',Interpreter='latex',FontSize=16);
    savefig('fig/figureLocNoAid.fig');
    print('-depsc2', 'fig/figureLocNoAid.eps');
    
    figure;scatter(stateNoSL(1).xGrid,stateNoSL(1).yGrid,'ks'); hold on;scatter(grid.x(indTarget),grid.y(indTarget),'r+');
    % % % scatter(sMGParam.grid.x(indTarget),sMGParam.grid.y(indTarget),'bs');
    legend('Grid location','Ground truth',Interpreter='latex',FontSize=15);xlim([-40,40]),ylim([-40,40]);xlabel('$x$ (m)',Interpreter='latex',FontSize=16);ylabel('$y$ (m)',Interpreter='latex',FontSize=16);
    savefig('fig/figureLocNoSL.fig');
    print('-depsc2', 'fig/figureLocNoSL.eps');
end

% scatter(stateProposed.gridRes.range,stateProposed.gridRes.azi);hold on;
%     gridGT = polar2cart(gridPolarGT);
% scatter(gridPolarGT.range,gridPolarGT.azi,'x');
%%
figure;
strLegend = strings(sum(baselineVec),1);
baselineInd = 1;

if baselineVec(3)
semilogy(1:2:Niter,mean(NMSENoSL(1:2:end,:),2),'->',LineWidth=1.5,Color = [138,43,226]/255);hold on;
strLegend(baselineInd) = 'OG-VAMP (w Fusion & Imperfect Tx prior)';
baselineInd = baselineInd+1;
end
if baselineVec(4)
semilogy(1:2:Niter,mean(NMSESV(1:2:end,:),2),'-<',LineWidth=1.5,Color = [255,128,0]/255);hold on;
strLegend(baselineInd) = 'JSSES (w/o Fusion)';
baselineInd = baselineInd+1;
end
if baselineVec(2)
semilogy(1:2:Niter,2*mean(NMSENoAid(1:2:end,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
strLegend(baselineInd) = 'JSSES (No aid)';baselineInd = baselineInd+1;
end
if baselineVec(2)
semilogy(1:2:Niter,mean(NMSEProposed(1:2:end,:),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
strLegend(baselineInd) = 'JSSES';
baselineInd = baselineInd+1;
end
if baselineVec(1)
semilogy(1:2:Niter,mean(NMSEOracle(1:2:end,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
strLegend(baselineInd) = 'Oracle';
baselineInd = baselineInd+1;
end
grid on;


% semilogy(NMSE,'-o');
% hold on;

% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend(strLegend,Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='best');
xlim([1 Niter]);
% ylim([3e-2,15]);
ylabel('NMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
% 
% nameStr = strcat('fig\figureConv','.fig');
% savefig(nameStr);
%%
figure;
strLegend = strings(sum(baselineVec),1);
baselineInd = 1;
if baselineVec(5)
semilogy(2:2:Niter,mean(RMSEVAMP(2:2:end,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
strLegend(baselineInd) = 'VAMP (w Imperfect Tx prior)';baselineInd = baselineInd+1;
end
if baselineVec(3)
semilogy(2:2:Niter,mean(RMSENoSL(2:2:end,:),2),'->',LineWidth=1.5,Color = [138,43,226]/255);hold on;
strLegend(baselineInd) = 'OG-VAMP (w Fusion & Imperfect Tx prior)';
baselineInd = baselineInd+1;
end
if baselineVec(4)
semilogy(2:2:Niter,mean(RMSESV(2:2:end,:),2),'-<',LineWidth=1.5,Color = [255,128,0]/255);hold on;
strLegend(baselineInd) = 'JSSES (w/o Fusion)';
baselineInd = baselineInd+1;
end
if baselineVec(2)
semilogy(2:2:Niter,mean(RMSENoAid(2:2:end,:),2),'-s',LineWidth=1.5,Color = [220,20,60]/255);hold on;
strLegend(baselineInd) = 'JSSES (No Aid)';
baselineInd = baselineInd+1;
semilogy(2:2:Niter,mean(RMSEProposed(2:2:end,:),2),'-o',LineWidth=1.5,Color = [0,0,255]/255);hold on;
strLegend(baselineInd) = 'JSSES';
baselineInd = baselineInd+1;
end
if baselineVec(1)
semilogy(2:2:Niter,mean(RMSEOracle(2:2:end,:),2),'-x',LineWidth=1.5,Color = [0,0,0]/255);hold on;
strLegend(baselineInd) = 'Oracle';
baselineInd = baselineInd+1;
end
grid on;


% semilogy(NMSE,'-o');
% hold on;

% legend('SBL','VBI','SVVBL','ALESD','OGMVSBLWOD','OOGMVSBL',Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='southwest');xlim([1 Niter]);

legend(strLegend,Interpreter='tex',FontName='Times New Roman',FontSize=10,Location='best');
xlim([1 Niter]);
ylim([1e-1,2]);
ylabel('RMSE',Interpreter='tex',FontName='Times New Roman',FontSize=15);
xlabel('Iteration index',Interpreter='tex',FontName='Times New Roman',FontSize=15);
if snr == 0
    nameStr = strcat('fig\figureConvRMSE','.fig');
    savefig(nameStr);
else
    nameStr = strcat('fig\figureConvRMSE2','.fig');
    savefig(nameStr);
end

%%

%%
str = datestr(now);
str = strrep(str,' ','');
str = strrep(str,'-','');
str = strrep(str,':','');

if snr == 0
str = strcat('results\',str,'mainOGConvtest','.mat');
save(str,"RMSENoAid",'RMSEProposed','RMSEOracle','RMSESV','RMSENoSL');
end
% save(str,"NMSENoAid",'NMSEProposed','NMSEOracle','NMSESV','NMSENoSL');
fprintf('\n Convergence Behavior Simulation Finished!\n');
%
