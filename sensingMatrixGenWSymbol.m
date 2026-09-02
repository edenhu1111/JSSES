function [Phi,outputDeriv] =...
    sensingMatrixGenWSymbol(xUAVinput,xR,subcarrInd,symbol,vecTO)
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
global Na Nv fc deltaF sigma2 c0 Nsize dArray biasAngle
Nr = size(xR,2);
Nt = size(xUAVinput,2);
% fn = fc + subcarrInd*deltaF;
fn = subcarrInd*deltaF;
Lambda = c0/fc;
fn = reshape(fn,[],1);  % reshape fn into a vector
Nsub = length(fn);
xGrid = linspace(-20,20,Nsize);
yGrid = linspace(-20,20,Nsize);
Phi = zeros(Nt*Nr*Na*Nv*Nsub,length(xGrid)*length(yGrid));
Phidx = zeros(size(Phi));
Phidy = zeros(size(Phi));
Phidz = zeros(size(Phi));

Phiddelay = zeros(Na*Nv*Nsub,length(xGrid)*length(yGrid),Nt*Nr);
PhiddelaywLoS = zeros(Nv*Na*Nsub,length(xGrid)*length(yGrid)+Nt*Nr,Nt*Nr);


for nt = 1:Nt
    xUAV = xUAVinput(:,nt);
    for ii = 1:length(xGrid)
        for jj = 1:length(yGrid)
            xTarget = [xGrid(ii);yGrid(jj);0];
            for nr = 1:Nr
                thetain = sin(asin((xTarget(2)-xR(2,nr))/norm(xTarget-xR(:,nr),2)) + biasAngle(nr));
                phiin = (xTarget(3)-xR(3,nr))/norm(xTarget-xR(:,nr),2);
                aa = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).');
                av = exp(-1j*2*pi/Lambda*dArray*phiin*(0:Nv-1).');
                aArray = kron(aa,av);
                tauT = norm(xUAV-xTarget,'fro')/c0;
                tauT2 = norm(xR(:,nr)-xTarget,'fro')/c0;

                sd = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr)));
                sdx = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr))).*...
                    (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(1)-xTarget(1));
                sdy = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr))).*...
                    (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(2)-xTarget(2));
                sdz = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr))).*...
                    (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(3)-xTarget(3));

                sdd = symbol.*exp(-1j*2*pi*fn*(tauT+tauT2+vecTO(nr))).*...
                    (-1j*2*pi*fn);


                rrtmp = kron(aArray,sd);
                rrtmpdx = kron(aArray,sdx);
                rrtmpdy = kron(aArray,sdy);
                rrtmpdz = kron(aArray,sdz);
                rrtmpddelay = kron(aArray,sdd);
                % Calculating Path Loss (eliminated)
                PL = 1;

                
                % Original Sensing Matrix Generation
                Phi((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    (ii-1)*length(yGrid) + jj) = PL*rrtmp;

                % The first-order Deriavative
                Phidx((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    (ii-1)*length(yGrid) + jj) = PL*rrtmpdx;
                Phidy((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    (ii-1)*length(yGrid) + jj) = PL*rrtmpdy;
                Phidz((nt-1)*Nr*Nv*Na*Nsub + (nr-1)*Nv*Na*Nsub + 1:(nt-1)*Nr*Nv*Na*Nsub + nr*Nv*Na*Nsub  ,...
                    (ii-1)*length(yGrid) + jj) = PL*rrtmpdz;
                Phiddelay(:,(ii-1)*length(yGrid) + jj,(nt-1)*Nr + nr) = PL*rrtmpddelay;
            end
        end
    end
end

%% LoS Component
% normPhi = vecnorm(Phi,2,1);
% normPhi = mean(normPhi);

LoSPhi = zeros(Nv*Na*Nsub,Nr,Nt);
LoSPhidx = zeros(Nv*Na*Nsub,Nr,Nt);
LoSPhidy = zeros(Nv*Na*Nsub,Nr,Nt);
LoSPhidz = zeros(Nv*Na*Nsub,Nr,Nt);
LoSPhiddelay = zeros(Nv*Na*Nsub,Nr,Nt);

for nt = 1:Nt
    xUAV = xUAVinput(:,nt);
    for nr = 1:Nr
        % rrtmp = []; rrtmpdx = [];rrtmpdy = [];rrtmpdz = [];rrtmpddelay = [];rrtmpdddd = [];

        thetain = sin(asin((xUAV(2)-xR(2,nr))/norm(xUAV-xR(:,nr),2)) + biasAngle(nr));
        phiin = (xUAV(3)-xR(3,nr))/norm(xUAV-xR(:,nr),2);
        aa = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).');
        av = exp(-1j*2*pi/Lambda*dArray*phiin*(0:Nv-1).');
        aArray = kron(aa,av);
        tauT = norm(xUAV-xR(:,nr),'fro')/c0;

        sd = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr)));
        sdx = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr))).*...
            (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(1)-xR(1,nr));
        sdy = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr))).*...
            (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(2)-xR(2,nr));
        sdz = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr))).*...
            (-1j*2*pi*fn)/(c0^2*tauT)*(xUAV(3)-xR(3,nr));

        sdd = symbol.*exp(-1j*2*pi*fn*(tauT+vecTO(nr))).*...
            (-1j*2*pi*fn);
        
        aadx = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).').*...
            (-1j*2*pi/Lambda*dArray*(0:Na-1).')*...
            (xUAV(1)-xR(1,nr))*(xUAV(2)-xR(2,nr))/norm(xTarget(1:2)-xR(1:2,nr),2)^3;
        aady = exp(-1j*2*pi/Lambda*dArray*thetain*(0:Na-1).').*...
            (-1j*2*pi/Lambda*dArray*(0:Na-1).')*...
            ((xUAV(1)-xR(1,nr))^2/norm(xTarget(1:2)-xR(1:2,nr),2)^3 + ...
            1/norm(xTarget(1:2)-xR(1:2,nr),2));

        rrtmp = kron(aArray,sd);
        rrtmpdx = kron(aArray,sdx)+kron(kron(aadx,av),sd);
        rrtmpdy = kron(aArray,sdy)+kron(kron(aady,av),sd);
        rrtmpdz = kron(aArray,sdz);
        rrtmpddelay = kron(aArray,sdd);
        % PL calculation
        PL = 1;
        rr = PL*rrtmp;
        nrm = 1;
        
        LoSPhi(:,nt,nr)   = rr/nrm;
        %1-st order derivative
        LoSPhidx(:,nt,nr) = PL*rrtmpdx/nrm;
        LoSPhidy(:,nt,nr) = PL*rrtmpdy/nrm;
        LoSPhidz(:,nt,nr) = PL*rrtmpdz/nrm;

        LoSPhiddelay(:,nt,nr) = PL*rrtmpddelay/nrm;


    end
end
LoSPhiDiag = [];LoSPhidxDiag=[];LoSPhidyDiag=[];LoSPhidzDiag=[];LoSPhiddelaywLoS=[];
for nt = 1:Nt
    for nr = 1:Nr
        LoSPhiDiag = blkdiag(LoSPhiDiag,LoSPhi(:,nt,nr)); 
        LoSPhidxDiag = blkdiag(LoSPhidxDiag,LoSPhidx(:,nt,nr)); 
        LoSPhidyDiag = blkdiag(LoSPhidyDiag,LoSPhidy(:,nt,nr)); 
        LoSPhidzDiag = blkdiag(LoSPhidzDiag,LoSPhidz(:,nt,nr)); 
        LoSPhiddelaywLoS = blkdiag(LoSPhiddelaywLoS,LoSPhiddelay(:,nt,nr));
    end
end
Phi = [LoSPhiDiag,Phi];
Phidx = [LoSPhidxDiag,Phidx];
Phidy = [LoSPhidyDiag,Phidy];
Phidz = [LoSPhidzDiag,Phidz];

for nt = 1:Nt
    for nr = 1:Nr
        PhiddelaywLoS(:,:,(nt-1)*Nr + nr) = ...
            [LoSPhiddelaywLoS((nt-1)*Nr*Na*Nv*Nsub + (nr-1)*Na*Nv*Nsub + 1:(nt-1)*Nr*Na*Nv*Nsub + nr*Na*Nv*Nsub,:),...
            Phiddelay(:,:,(nt-1)*Nr + nr)];
    end
end
nrm = norm(Phi,'fro');
Phi = Phi/nrm;
outputDeriv.Phidx = Phidx/nrm;
outputDeriv.Phidy = Phidy/nrm;
outputDeriv.Phidz = Phidz/nrm;
outputDeriv.PhiddelaywLoS = PhiddelaywLoS/nrm;


end