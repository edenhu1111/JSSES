 function [Phi,outputDeriv] =...
    sensingMatrixGenWSymbolPolar(sMGParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: sensingMatrixGen
% Description: Generate sensing matrix for compressed sensing based imaging
% Author: Eden Hu (huyb@mail.sim.ac.cn)
%
% Input description
% xUav: UAV Position(Estimated)
% xR: BS Position(Accurate)
% subcarrInd: used subcarriers
% vecTO: time offset vector(Estimated)
%
% Output description
% Phi: generated sensing matrix
%               Nr->Na->Nv->Nsub
%               2D->1D: x outer, y inner (similar to vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flagLoS = sMGParam.flagLoS;

xUAVinput = sMGParam.xUAVinput;
xR= sMGParam.xR;
subcarrInd = sMGParam.subcarrInd;
symbol = sMGParam.symbol;
vecTO = sMGParam.vecTO;
grid = sMGParam.grid;  % Polar

Na = sMGParam.Na;
Nv = sMGParam.Nv;
fc = sMGParam.fc;
deltaF = sMGParam.deltaF;
c0 = physconst('LightSpeed');
dArray = sMGParam.dArray;
ea = sMGParam.ULALine;
ev = [0;0;1];
% global Na Nv fc deltaF c0 dArray biasAngle
Nr = size(xR,2);
Nt = size(xUAVinput,2);
% fn = fc + subcarrInd*deltaF;
fn = subcarrInd*deltaF;
Lambda = c0/fc;
fn = reshape(fn,[],1);  % reshape fn into a vector
Nsub = length(fn);
% xGrid = linspace(-20,20,Nsize);
% yGrid = linspace(-20,20,Nsize);
xGrid = grid.range.*cos(grid.azi);
yGrid = grid.range.*sin(grid.azi);
% xGrid = grid.x;
% yGrid = grid.y;

Phi = zeros(Nt*Nr*Na*Nv*Nsub,length(xGrid));
Phidx = zeros(size(Phi,1),size(Phi,2));
Phidy = zeros(size(Phi,1),size(Phi,2));
if flagLoS
    PhidxUAV = zeros(size(Phi,1),size(Phi,2)+1);
    PhidyUAV = zeros(size(Phi,1),size(Phi,2)+1);
    PhidzUAV = zeros(size(Phi,1),size(Phi,2)+1);
    Phidtau  = zeros(size(Phi,1),size(Phi,2)+1);
else
    PhidxUAV = zeros(size(Phi,1),size(Phi,2));
    PhidyUAV = zeros(size(Phi,1),size(Phi,2));
    PhidzUAV = zeros(size(Phi,1),size(Phi,2));
    Phidtau  = zeros(size(Phi,1),size(Phi,2));
end


for nt = 1:Nt
    xUAV = xUAVinput(:,nt);
    for ii = 1:length(xGrid)
            xTarget = [xGrid(ii);yGrid(ii);0];
            for nr = 1:Nr
                thetain = (xTarget-xR(:,nr))'*ea(:,nt,nr)/norm(xTarget-xR(:,nr),'fro');
                phiin =   (xTarget-xR(:,nr))'*ev/norm(xTarget-xR(:,nr),'fro');
                aa = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).');
                av = exp(-1j*2*pi/Lambda*dArray*phiin*(0:Nv-1).');
                aArray = kron(aa,av);

                dT1 = norm(xUAV-xTarget,'fro');
                dT2 = norm(xR(:,nr)-xTarget,'fro');
                % tauT = norm(xUAV-xTarget,'fro')/c0;
                tauT  = dT1/c0;
                tauT2 = dT2/c0;

                aadx = (1j*pi*(0:Na-1)').*((xTarget-xR(:,nr))'*ea(:,nt,nr).*(xTarget(1)-xR(1,nr)) - ea(1,nt,nr)*dT2^2)/dT2^3.*aa;
                aady = (1j*pi*(0:Na-1)').*((xTarget-xR(:,nr))'*ea(:,nt,nr).*(xTarget(2)-xR(2,nr)) - ea(2,nt,nr)*dT2^2)/dT2^3.*aa;

                avdx = (1j*pi*(0:Nv-1)').*((xTarget-xR(:,nr))'*ev.*(xTarget(1)-xR(1,nr)) - ev(1)*dT2^2)/dT2^3.*av;
                avdy = (1j*pi*(0:Nv-1)').*((xTarget-xR(:,nr))'*ev.*(xTarget(2)-xR(2,nr)) - ev(2)*dT2^2)/dT2^3.*av;

                sd = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr)));
                sdx = (-1j*2*pi*fn/c0)*( (xTarget(1)-xUAV(1))/dT1 + (xTarget(1)-xR(1,nr))/dT2 ).*sd;
                sdy = (-1j*2*pi*fn/c0)*( (xTarget(2)-xUAV(2))/dT1 + (xTarget(2)-xR(2,nr))/dT2 ).*sd;
                
                sdxU = (-1j*2*pi*fn/c0)*( (xUAV(1) - xTarget(1))/dT1).*sd;
                sdyU = (-1j*2*pi*fn/c0)*( (xUAV(2) - xTarget(2))/dT1).*sd;
                sdzU = (-1j*2*pi*fn/c0)*( (xUAV(3) - xTarget(3))/dT1).*sd;

                sdtau = (-1j*2*pi*fn).*sd;



                rrtmp = kron(aArray,sd);
                rrtmpdx =   kron(kron(aadx,av),sd) + kron(kron(aa,avdx),sd) + kron(aArray,sdx);
                rrtmpdy =   kron(kron(aady,av),sd) + kron(kron(aa,avdy),sd) + kron(aArray,sdy);

                rrtmpdxU = kron(aArray,sdxU);
                rrtmpdyU = kron(aArray,sdyU);
                rrtmpdzU = kron(aArray,sdzU);

                rrtmpdtau = kron(aArray,sdtau);


                % Calculating Path Loss (eliminated)
                PL = 1;

                
                % Original Sensing Matrix Generation
                Phi((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii) = PL*rrtmp;

                % The first-order Deriavative
                Phidx((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdx;
                Phidy((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdy;
                PhidxUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdxU;
                PhidyUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdyU;
                PhidzUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdzU;
                Phidtau((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    ii+flagLoS) = PL*rrtmpdtau;
            end
    end
end


%% LoS Component
% normPhi = vecnorm(Phi,2,1);
% normPhi = mean(normPhi);
if flagLoS
LoSPhi = zeros(Nv*Na*Nsub,Nr,Nt);

for nt = 1:Nt
    xUAV = xUAVinput(:,nt);
    for nr = 1:Nr
        % rrtmp = []; rrtmpdx = [];rrtmpdy = [];rrtmpdz = [];rrtmpddelay = [];rrtmpdddd = [];

        thetain = ((xUAV-xR(:,nr))'*ea(:,nt,nr))/(norm(xUAV-xR(:,nr),2)*norm(ea(:,nt,nr),'fro'));
        phiin =   ((xUAV-xR(:,nr))'*ev)/norm(xUAV-xR(:,nr),2);
        aa = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).');
        av = exp(-1j*2*pi/Lambda*dArray*phiin*(0:Nv-1).');
        aArray = kron(aa,av);
        tauT = norm(xUAV-xR(:,nr),'fro')/c0;

        sd = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr)));

        sdxU = (-1j*2*pi*fn/c0)*( (xUAV(1) - xR(1,nr))/dT1).*sd;
        sdyU = (-1j*2*pi*fn/c0)*( (xUAV(2) - xR(2,nr))/dT1).*sd;
        sdzU = (-1j*2*pi*fn/c0)*( (xUAV(3) - xR(3,nr))/dT1).*sd;

        sdtau = (-1j*2*pi*fn).*sd;

        rrtmp = kron(aArray,sd);

        rrtmpdxU = kron(aArray,sdxU);
        rrtmpdyU = kron(aArray,sdyU);
        rrtmpdzU = kron(aArray,sdzU);

        rrtmpdtau = kron(aArray,sdtau);
        % PL calculation
        PL = 1;
        rr = PL*rrtmp;
        nrm = 1;
        
        LoSPhi(:,nt,nr)   = rr/nrm;
        PhidxUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
            1) = PL*rrtmpdxU;
        PhidyUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
            1) = PL*rrtmpdyU;
        PhidzUAV((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
            1) = PL*rrtmpdzU;
        Phidtau((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
            1) = PL*rrtmpdtau;

    end
end
LoSPhiDiag = [];
for nt = 1:Nt
    for nr = 1:Nr
        LoSPhiDiag = blkdiag(LoSPhiDiag,LoSPhi(:,nt,nr)); 
    end
end
Phi = [LoSPhiDiag,Phi];
end

Phidr   =  Phidx.*repmat(cos(grid.azi),Nt*Nr*Na*Nv*Nsub,1) + Phidy.*repmat(sin(grid.azi),Nt*Nr*Na*Nv*Nsub,1);
Phidazi = (-Phidx.*repmat(sin(grid.azi),Nt*Nr*Na*Nv*Nsub,1) +...
    Phidy.*repmat(cos(grid.azi),Nt*Nr*Na*Nv*Nsub,1)).*repmat(grid.range,Nt*Nr*Na*Nv*Nsub,1);


nrm = norm(Phi,'fro');
Phi = Phi/nrm;
% outputDeriv.Phidx = Phidx/nrm;
% outputDeriv.Phidy = Phidy/nrm;
outputDeriv.Phidr = Phidr/nrm;
outputDeriv.Phidazi = Phidazi/nrm;

outputDeriv.PhidxUAV = PhidxUAV/nrm;
outputDeriv.PhidyUAV = PhidyUAV/nrm;
outputDeriv.PhidzUAV = PhidzUAV/nrm;
outputDeriv.Phidtau  = Phidtau/nrm;

end