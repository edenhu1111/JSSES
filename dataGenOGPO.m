%% ISAC Tx & Rx
clear all;
global Nsize 
fprintf('Data Generator\n');

load('dataset/xOG');
load('dataset/1/grid');
load('dataset/1/xUAV');

load('dataset\1\xR','xR');
load('dataset\1\ULALine','ULALine');
% System parameters
c0 = physconst('lightspeed');  % Speed of light
fc = 28e+9;                    % Carrier frequency (Ka band)
% BW = 60e6;
Ns = 1024;
deltaF = 30e3;   % Subcarrier spacing
BW = Ns*deltaF;
% BW = deltaF * Ns;  % Total bandwidth
T = 1 / deltaF;    % Symbol duration
Tcp = T / 4;       % Cyclic prefix duration
Ts = T + Tcp;      % Total symbol duration including cyclic prefix
sigmaPos = 1/sqrt(3);
sigma2 = 10^(-60/10);
rRange = 80;       % Range limit of interest in meters
aziRange = pi/2;
niter = 20;        % Number of iterations for estimation
Na = 32;           % Number of antenna elements
Nv = 1;            % Number of virtual elements (for array processing)
Nsize = 20;        % Size parameter (likely for some internal buffer or iteration limit)
NAzi = 16;
NRange = 16;

Nsim = 100;        % Number of simulation runs

% Receiver position [x; y; z]
biasY = 10;
% xR = [  -rRange/2,         0,     -3*rRange/8,  -3*rRange/8,        0 ,      -rRange/2,    -rRange/2;...
%             0,         -rRange/2,           +rRange/2,    -rRange/2,      rRange/2,   +rRange/2,    -rRange/2;
%             0,             0,              0,             0,             0 ,          0  ,           0    ];  % Example receiver position at x = -50, y = 0, z = 5
% ULALine = [0, 1, 1, 1, 1, 1;  ...
%            1, 0, 0, 0, 1, -1; 
%            0, 0 ,0, 0, 0, 0];
% xUAVGT = [rRange/2;rRange/2;0];
xUAVGT = xUAV;

xR = xR(:,[1:2]);
ULALine = ULALine(:,[1:2]);




Nr = size(xR,2);

% Number of receivers

% Pilot and data subcarrier indices
pilotSubcarrInd = 0:32:Ns-1;
% dataSubcarrInd = [];
% scaleFac = 5;
% for ii = 1:scaleFac
%     dataSubcarrInd = [dataSubcarrInd , pilotSubcarrInd + ii];
% end
% dataSubcarrInd = sort(dataSubcarrInd);
% dataNum = length(dataSubcarrInd);
arrayPO = [2,4,6,8,10];
xing = [1,0,-1, 0;
        0,1, 0,-1;
        0,0, 0, 0];
for indPO = 1:length(arrayPO)
    for indXing = 1:4
        xUAVGT(:,indPO,indXing) = xUAV + xing(:,indXing)*arrayPO(indPO);
    end
end
pilot = ones(length(pilotSubcarrInd),1);

flagBias = 1;
% xUAV = xUAVGT + flagBias*(2*rand(size(xUAVGT))-1);
% xUAV(3) = xUAVGT(3);
tauGT = 0.25/BW;
% for ii = 1:length(tauGT)
%     tauGT(ii) = tauGT(ii)/BW*flagBias;
% end
tauEST = 0;

dArray = c0/fc/2;
biasAngle = [0]';                                                         %degree

% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
xTarget = zeros(3,9);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation
PhipGeniePO = cell(length(tauGT),1);
PhiGPO = cell(length(tauGT),1);
derivGPO = cell(length(tauGT),1);
sMGParamGPO = cell(length(tauGT),1);
zpPO = cell(length(tauGT),1);

for ii = 1:length(arrayPO)
    for jj = 1:4
        Phip = zeros(Na*Nv*length(pilotSubcarrInd),Nsize^2,Nr);
        
        PhiG = Phip;
        PhipUni = Phip;
        % PhidUni = Phid;
        
        yp = zeros(size(Phip,1),Nr);
        
        
        
        % numTarget = 10;
        indTmp = [(Nsize)/4+1,2*(Nsize)/4+1,3*(Nsize)/4+1];
        % indTmp = [(NAzi-1)/4+1,(NAzi-1)/4+1,3*(NAzi-1)/4+1];
        
        indTmp = round(indTmp);
        targetInd = [indTmp + Nsize*(round((Nsize)/4)),indTmp + 2*Nsize*(round((Nsize)/4)),indTmp + 3*Nsize*(round((Nsize)/4))];
        targetInd = targetInd([1:2,4:5,9]);
        
        
        numTarget = numel(targetInd);
        
        
        % x = zeros(Nsize^2,Nr);
        
        eyeS = eye(Na*Nv);
        
        deltaGrid = (rRange)/(Nsize-1);
        % 
        % range = linspace(40,rRange,NRange);
        % azi   = linspace(-aziRange/2,aziRange/2,NAzi);
        % gridPolar.azi = kron(ones(size(range)),azi);
        % gridPolar.range = kron(range,ones(size(azi)));
        % 
        % grid.x = gridPolar.range.*cos(gridPolar.azi) - 40;
        % grid.y = gridPolar.range.*sin(gridPolar.azi);
        
        indGrid = linspace(-rRange/2,rRange/2,Nsize+2); indGrid = indGrid(2:end-1);
        
        gridUni.x = kron(indGrid,ones(1,Nsize));
        gridUni.y =  kron(ones(1,Nsize),indGrid);
        
        
        sMGParam.xR = xR;
        % sMGParam.grid = grid;
        sMGParam.Na = Na;
        sMGParam.Nv = Nv;
        sMGParam.fc = fc;
        sMGParam.deltaF = deltaF;
        sMGParam.dArray = dArray;
        sMGParam.biasAngle = biasAngle;
        sMGParam.flagLoS = 0;
        sMGParam.subcarrInd = pilotSubcarrInd;
        sMGParam.NRange = NRange;sMGParam.NAzi = NAzi;
        
        
        tic
        for tt = 1:Nr
        %     PhipUni(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,pilotSubcarrInd,pilot,  tauGT);
        %     PhidUni(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,symbols,tauGT);
            
        %     PhidUni(:,:,tt) = sensingMatrixGenWSymbolOG(sMGParam);
            sMGParam.xUAVinput = xUAVGT(:,ii,jj);
            sMGParam.xR = xR(:,tt);
            sMGParam.symbol = pilot;
            sMGParam.ULALine = ULALine(:,tt);
        %     [Phip(:,:,tt),derivP] = sensingMatrixGenWSymbolOG(xUAV(:,tt),xR,pilotSubcarrInd,pilot,  tauGT,grid);
            sMGParam.grid = grid;
            sMGParam.vecTO = tauGT;
            [Phip(:,:,tt),~] = sensingMatrixGenWSymbolOG(sMGParam);
        
            sMGParamG = sMGParam;
            sMGParamG.grid = gridUni;
            sMGParamG.xUAVinput = xUAVGT(:,ii,jj);
            sMGParamG.vecTO = tauGT;
            
            [PhiG(:,:,tt),derivG(tt)]   = sensingMatrixGenWSymbolOG(sMGParamG);
            if ii == 1 && jj == 1
                sMGParam.grid = gridUni;
                sMGParam.xUAVinput = xUAV;
                sMGParam.vecTO = tauEST;
            
                [PhipUni(:,:,tt),derivP(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
            end
        
        %     if ii == 1
        %         PhipUni(:,:,tt) = PhipUni(:,:,tt)*sqrt(size(PhipUni,2));
        %     
        %         derivP(tt).Phidx = derivP(tt).Phidx*sqrt(size(PhipUni,2));
        %         derivP(tt).Phidy = derivP(tt).Phidy*sqrt(size(PhipUni,2));
        %     %     derivP(tt).Phidr   = derivP(tt).Phidr*sqrt(size(PhipUni,2));
        %     %     derivP(tt).Phidazi = derivP(tt).Phidazi*sqrt(size(PhipUni,2));
        %         derivP(tt).PhidxUAV = derivP(tt).PhidxUAV*sqrt(size(PhipUni,2));
        %         derivP(tt).PhidyUAV = derivP(tt).PhidyUAV*sqrt(size(PhipUni,2));
        %         derivP(tt).PhidzUAV = derivP(tt).PhidzUAV*sqrt(size(PhipUni,2));
        %         derivP(tt).Phidtau  = derivP(tt).Phidtau*sqrt(size(PhipUni,2));
        %     end
        %     Phip(:,:,tt) = Phip(:,:,tt)*sqrt(size(Phip,2));
        %     PhiG(:,:,tt) = PhiG(:,:,tt)*sqrt(size(PhiG,2));
        %     derivG(tt).Phidx = derivG(tt).Phidx*sqrt(size(PhipUni,2));
        %     derivG(tt).Phidy = derivG(tt).Phidy*sqrt(size(PhipUni,2));
        % %     derivG(tt).Phidr = derivG(tt).Phidr*sqrt(size(PhipUni,2));
        % %     derivG(tt).Phidazi = derivG(tt).Phidazi*sqrt(size(PhipUni,2));
        %     derivG(tt).PhidxUAV = derivG(tt).PhidxUAV*sqrt(size(PhipUni,2));
        %     derivG(tt).PhidyUAV = derivG(tt).PhidyUAV*sqrt(size(PhipUni,2));
        %     derivG(tt).PhidzUAV = derivG(tt).PhidzUAV*sqrt(size(PhipUni,2));
        %     derivG(tt).Phidtau  = derivG(tt).Phidtau*sqrt(size(PhipUni,2));
        %     
        
            zp(:,tt) = Phip(:,:,tt)*x(:,tt);
        end
        if ii == 1 && jj == 1
            sMGParam.grid = gridUni;
            sMGParam.derivP = derivP;
        end
        sMGParamG.derivP = derivG;
        sMGParamG.grid = gridUni;
          
        PhipGeniePO{ii,jj} = Phip;
        PhiGPO{ii,jj} = PhiG;
        derivGPO{ii,jj} = derivG;
        sMGParamGPO{ii,jj} = sMGParamG;
        
        zpPO{ii,jj} = zp;
    end
end
toc






save('dataset\3\PhipGenie','PhipGeniePO');
% save('dataset\1\Phip','Phip');
save('dataset\3\PhiG','PhiGPO');


save('dataset\3\tauGT',"tauGT");
% save('dataset\1\xR','xR');
% save('dataset\1\ULALine','ULALine');

save('dataset\3\zp','zpPO');


% save('dataset\2\derivP','derivP');
save('dataset\3\derivG','derivGPO');
% save('dataset\2\sMGParam','sMGParam');
save('dataset\3\sMGParamG','sMGParamGPO');

% save('dataset\1\grid',"grid");