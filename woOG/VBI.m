function [xest, state] = VBI(yp,Phip,inputParam)
%%   Multi-view variation bayesian inference
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
tic
global gg
%% Initialization
a    = 1e0;  b = 1e0;
abar = 1e0;   bbar = 1e-8;
c = 1e-8;   d = 1e-8;

Mp = size(Phip,1);
N = size(Phip,2);
T = size(Phip,3);
eps = 1e-6;
convBreak = inputParam.convBreak;
gammaOmega = inputParam.gammaOmega*ones(T,1);
Niter = inputParam.Niter;
NiterA = inputParam.NiterA;
lambda0 = inputParam.lambda0;
% xUAV = inputParam.xUAV;
% xR = inputParam.xR;
% dataSubcarrInd = inputParam.dataSubcarrInd;
% pilotSubcarrInd = inputParam.pilotSubcarrInd;
% lambdaIn = 0.5*ones(N-1,T);


% d = zeros(dataNum,T);
st = lambda0*ones(N,T);
st(1,:) = ones(1,T);

rho = (a.*st + abar.*(1-st))./(b.*st + bbar.*(1-st)).*ones(N,T);
xest = zeros(N,T);


lambdaIn = lambda0*ones(N-1,T);

% lnrho = zeros(N,T);
nmse = zeros(NiterA,Niter);
%% VBI Iteration
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmega(tt)...
                 + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt) )*gammaOmega(tt);

            muz = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            ss = abs(diag(Sx(:,:,tt)));    

            %% update rho and gammaOmega
            rho(:,tt) =  (st(:,tt).*a + (1-st(:,tt)).*abar + 1)./...
                (abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
%             rho(:,tt) =  (a + 1)./...
%                 (abs(xest(:,tt)).^2 + ss + b); 
            gammaOmega(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muz +d);

            lnrho(:,tt) = psi(st(:,tt).*a + (1-st(:,tt)).*abar + 1) - log(abs(xest(:,tt)).^2 + ss + st(:,tt).*b + (1-st(:,tt)).*bbar);
            %% update st
    
%             st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt))./...
%                 (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)) + ...
%                 bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)));
            st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
                (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
                bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
        end
        %% Break when converging
        % to be completed
        nmse(niterA,niter) = mean(vecnorm(xest-gg,2).^2./vecnorm(gg,2).^2);
        if norm(xest - xold,'fro')^2/norm(xold,'fro')^2 < 1e-5 && convBreak
            break;
        end
        %
        %%%%%%%%%%%%%%%%%
    end
    % Break when converging
    if norm(xest - xoldOut,'fro')^2/norm(xold,'fro')^2 < 1e-5 && convBreak
        break;
    end
end
state.NMSE = nmse;
state.rho = rho;
state.gammaOmega = gammaOmega;
state.telap = toc;
end