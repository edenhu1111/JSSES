%% ISAC Tx & Rx
clear all;
global Nsize 
fprintf('Data Generator\n');
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
xR = [  -rRange/2,         0,     -3*rRange/8,  -3*rRange/8,        0 ,      -rRange/2,    -rRange/2;...
            0,         -rRange/2,           +rRange/2,    -rRange/2,      rRange/2,   +rRange/2,    -rRange/2;
            0,             0,              0,             0,             0 ,          0  ,           0    ];  % Example receiver position at x = -50, y = 0, z = 5
ULALine = [0, 1, 1, 1, 1, 1;  ...
           1, 0, 0, 0, 1, -1; 
           0, 0 ,0, 0, 0, 0];
% xUAVGT = [rRange/2;rRange/2;0];
xUAVGT = [-rRange/2-5;-rRange/2-5;0];

xR = xR(:,[1:2]);
ULALine = ULALine(:,[1:2]);
% xR = [-rRange/2,rRange/2,      0,             0 ;...
%         0,          0 ,     -rRange/2,  rRange/2;
%         5,          5,         5,             5];  % Example receiver position at x = -50, y = 0, z = 5
% ULALine = [0, 0, 1, 1;  ...
%            1, 1 ,0, 0; 
%            0, 0 ,0, 0];



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

pilot = ones(length(pilotSubcarrInd),1);

flagBias = 1;
xUAV = xUAVGT + flagBias*(2*rand(size(xUAVGT))-1);
xUAV(3) = xUAVGT(3);

tauGT = 0.25/BW*flagBias;
tauEST = zeros(size(tauGT));

dArray = c0/fc/2;
biasAngle = [0]';                                                         %degree

% xTarget = [-50,0,50,-50,0,50,-50,0,50; -50,-50,-50,0,0,0,50,50,50; 0,0,0,0,0,0,0,0,0];
xTarget = zeros(3,9);
% xTarget = [5,5,-5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;5,-5,5,-5,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5;0,0,0,0,5,5,-5,-5,5,5,-5,-5,5,5,-5,-5];
%% Signal generation
% Phip = zeros(Na*Nv*length(pilotSubcarrInd),NAzi*NRange,Nr);
Phip = zeros(Na*Nv*length(pilotSubcarrInd),Nsize^2,Nr);

PhiG = Phip;
PhipUni = Phip;
% PhidUni = Phid;

yp = zeros(size(Phip,1),Nr);
% yd = zeros(size(Phid,1),Nr);



% numTarget = 10;
indTmp = [(Nsize)/4+1,2*(Nsize)/4+1,3*(Nsize)/4+1];
% indTmp = [(NAzi-1)/4+1,(NAzi-1)/4+1,3*(NAzi-1)/4+1];

indTmp = round(indTmp);
% targetInd = [indTmp + NAzi*round((NRange-1)/4),indTmp + 2*NAzi*round((NRange-1)/4),indTmp + 3*NAzi*round((NRange-1)/4)];
targetInd = [indTmp + Nsize*(round((Nsize)/4)),indTmp + 2*Nsize*(round((Nsize)/4)),indTmp + 3*Nsize*(round((Nsize)/4))];
targetInd = targetInd([1:2,4:5,9]);
% targetInd = targetInd(5);
% targetInd = 212;
% targetInd = round(targetInd);
numTarget = numel(targetInd);
% targetInd = (Nsize-1)/4+1

% PhidwithSym = Phid;
% d = zeros(dataNum,modNum);
% x = zeros(NRange*NAzi,Nr);
x = zeros(Nsize^2,Nr);

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
grid.x = kron(indGrid,ones(1,Nsize));
grid.y = kron(ones(1,Nsize),indGrid);
gridUni.x = grid.x;
gridUni.y = grid.y;

% gridPolarGT = gridPolar;
% gridPolarGT.azi(targetInd) = gridPolar.azi(targetInd) + 0.8*(rand(1,length(targetInd))-0.5)*2*pi/3/(NAzi-1);
% gridPolarGT.range(targetInd) = gridPolar.range(targetInd) + 0.8*(rand(1,length(targetInd))-0.5)*rRange/(NRange-1);


% grid.x = gridPolarGT.range.*cos(gridPolarGT.azi) - 40;
% grid.y = gridPolarGT.range.*sin(gridPolarGT.azi);

grid.x(targetInd) = grid.x(targetInd) + 0.8*deltaGrid * (2*rand(1,length(targetInd))-1)/2;
grid.y(targetInd) = grid.y(targetInd) + 0.8*deltaGrid * (2*rand(1,length(targetInd))-1)/2;

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
    sMGParam.xUAVinput = xUAVGT;
    sMGParam.xR = xR(:,tt);
    sMGParam.symbol = pilot;
    sMGParam.ULALine = ULALine(:,tt);
%     [Phip(:,:,tt),derivP] = sensingMatrixGenWSymbolOG(xUAV(:,tt),xR,pilotSubcarrInd,pilot,  tauGT,grid);
    sMGParam.grid = grid;
    sMGParam.vecTO = tauGT;
    [Phip(:,:,tt),~] = sensingMatrixGenWSymbolOG(sMGParam);

    sMGParamG = sMGParam;
    sMGParamG.grid = gridUni;
    sMGParamG.xUAVinput = xUAVGT;
    sMGParamG.vecTO = tauGT;
    
    [PhiG(:,:,tt),derivG(tt)]   = sensingMatrixGenWSymbolOG(sMGParamG);

    sMGParam.grid = gridUni;
    sMGParam.xUAVinput = xUAV;
    sMGParam.vecTO = tauEST;

    [PhipUni(:,:,tt),derivP(tt)] = sensingMatrixGenWSymbolOG(sMGParam);


%     sMGParam.subcarrInd = dataSubcarrInd;
%     sMGParam.symbol = symbols;

%     sMGParam.grid = grid;
%     sMGParam.xUAVinput = xUAVGT;
%     sMGParam.vecTO = tauGT;
%     
% %     [Phid(:,:,tt),~] = sensingMatrixGenWSymbolOG(sMGParam);
% %     Phid(:,:,tt) = Phid(:,2:end,tt);
%     sMGParam.grid = gridUni;
%     sMGParam.xUAVinput = xUAV;
%     sMGParam.vecTO = tauEST;
%     [PhidUni(:,:,tt),derivD(tt)] = sensingMatrixGenWSymbolOG(sMGParam);
%     [Phid(:,:,tt),derivD] = sensingMatrixGenWSymbolOG(xUAV(:,tt),xR,dataSubcarrInd ,symbols,tauGT,grid);

%     Phid(:,:,tt) = Phid(:,:,tt)*sqrt(size(Phid,1))/sqrt(size(Phip,1));
%     PhipUni(:,:,tt) = PhipUni(:,:,tt)*sqrt(size(PhipUni,2));
%     Phip(:,:,tt) = Phip(:,:,tt)*sqrt(size(Phip,2));
%     PhiG(:,:,tt) = PhiG(:,:,tt)*sqrt(size(PhiG,2));
%     derivP(tt).Phidx = derivP(tt).Phidx*sqrt(size(PhipUni,2));
%     derivP(tt).Phidy = derivP(tt).Phidy*sqrt(size(PhipUni,2));
% 
%     derivP(tt).PhidxUAV = derivP(tt).PhidxUAV*sqrt(size(PhipUni,2));
%     derivP(tt).PhidyUAV = derivP(tt).PhidyUAV*sqrt(size(PhipUni,2));
%     derivP(tt).PhidzUAV = derivP(tt).PhidzUAV*sqrt(size(PhipUni,2));
%     derivP(tt).Phidtau  = derivP(tt).Phidtau*sqrt(size(PhipUni,2));
% 
%     derivG(tt).Phidx = derivG(tt).Phidx*sqrt(size(PhipUni,2));
%     derivG(tt).Phidy = derivG(tt).Phidy*sqrt(size(PhipUni,2));
%     derivG(tt).PhidxUAV = derivG(tt).PhidxUAV*sqrt(size(PhipUni,2));
%     derivG(tt).PhidyUAV = derivG(tt).PhidyUAV*sqrt(size(PhipUni,2));
%     derivG(tt).PhidzUAV = derivG(tt).PhidzUAV*sqrt(size(PhipUni,2));
%     derivG(tt).Phidtau  = derivG(tt).Phidtau*sqrt(size(PhipUni,2));

    x(targetInd,tt) =  (5+randn(numTarget,1)).*exp(1j*2*pi*rand(numTarget,1));

%     x(1,tt) = 8 + 10/sqrt(2)*(randn + 1j*randn);
%     x(1,tt) = 1/sqrt(2)*(randn + 1j*randn);
%     x(1,tt) = 0*exp(1j*2*pi*rand);

    zp(:,tt) = Phip(:,:,tt)*x(:,tt);
%     for nMod = 1:modNum
%         PhidwithSym(:,:,tt,nMod) = repmat(kron(ones(Na*Nv,1),d(:,tt,nMod)),1,size(Phid,2)).*Phid(:,:,tt);
%         zd(:,tt,nMod) = PhidwithSym(:,:,tt,nMod)*x(:,tt);
%     end
end
toc

noisedatap = (randn(size(Phip,1),Nr,Nsim) + 1j*randn(size(Phip,1),Nr,Nsim))/sqrt(2);
% noisedatad = (randn(size(Phid,1),Nr,Nsim) + 1j*randn(size(Phid,1),Nr,Nsim))/sqrt(2);
gg = x;

PhipGenie = Phip;
% PhidGenie = Phid;

Phip = PhipUni;
% Phid = PhidUni;
sMGParam.grid = gridUni;
sMGParam .derivP = derivP;
sMGParamG.grid = gridUni;
sMGParamG.derivP = derivG;
save('dataset\1\noisedatap','noisedatap');
% save('dataset\1\noisedatad','noisedatad');
save('dataset\1\PhipGenie','PhipGenie');
% save('dataset\1\PhidGenie','PhidGenie');
% save('dataset\1\Phid','Phid');
save('dataset\1\Phip','Phip');
save('dataset\1\PhiG','PhiG');

save('dataset\1\xUAVGT','xUAVGT');
save('dataset\1\xUAV','xUAV');
save('dataset\1\tauGT',"tauGT");
save('dataset\1\xR','xR');
save('dataset\1\ULALine','ULALine');

save('dataset\1\zp','zp');
% save('dataset\1\zd','zd');
% save('dataset\1\d','d');

save('dataset\xOG','x');
save('dataset\1\derivP','derivP');
save('dataset\1\derivG','derivG');
save('dataset\1\sMGParam','sMGParam');
save('dataset\1\sMGParamG','sMGParamG');

save('dataset\1\grid',"grid");