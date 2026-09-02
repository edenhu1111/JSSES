function [xest,dPost, state] = VBI2(yp,yd,Phip,Phid, inputParam)
%%   Multi-view variation bayesian inference
%   y: measurements, Phi: measurement matrix, sigma: noise std. The
%   algorithm is executed till ||y - Phi x|| < eps if Nmax = 0 else Nmax iterations. Returns sparse x and
%   time elapsed
% tic
global gg dd
% gammaDF = @(x,a,b)(x.^(a-1).*exp(-b.*x).*b^(a)./gamma(a));

%% Initialization
constell = inputParam.constell;
% a    = 1e0;  b = 1e1;
% abar = 1e0;   bbar = 1e-8;

a = 1;   b = 1e-3; 
c = 1e-3;   d = 1e-3;
% abar = 1e0;   bbar = 1e-6;
% if SV
%     abar = 1e0;   bbar = 1e-6;
% else
% %     abar = 1e-50;   bbar = 1e-50;
% end


Mp = size(Phip,1);
Md = size(Phid,1);
N = size(Phip,2);
T = size(Phip,3);
eps = 1e-3;
convBreak = inputParam.convBreak;
gammaOmega = inputParam.gammaOmega*ones(T,1);
gammaOmegaD = gammaOmega;
gammaOmegaP = gammaOmega;
Niter = inputParam.Niter;
NiterA = inputParam.NiterA;
dataNum = inputParam.dataNum;
JESD = inputParam.JESD;
lambda0 = inputParam.lambdaS*inputParam.lambda0;

% xUAV = inputParam.xUAV;
% xR = inputParam.xR;
% dataSubcarrInd = inputParam.dataSubcarrInd;
% pilotSubcarrInd = inputParam.pilotSubcarrInd;
% lambdaIn = 0.5*ones(N-1,T);

NaNv = size(inputParam.eyeS,1);

% d = zeros(dataNum,T);

rho = 1e2*ones(N,T);
% rho = ones(N,T);

% rho = (a)./(b).*ones(N,T);

xest = zeros(N,T);
PhidPost = zeros(size(Phid));
% PhidPostSquared = zeros(size(Phid(:,:,1)'*Phid(:,:,1)));
PhidPostSquared = zeros(size(Phid,2),size(Phid,2),T);
v  = zeros(dataNum,T);
mu = zeros(dataNum,T);
dPost = zeros(dataNum,T);
sigmaDPost = zeros(dataNum,T);

% lnrho = zeros(N,T);
nmse = zeros(NiterA,Niter);

if 1
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T  
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                 + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt));
            muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
            rho(:,tt) = (a + 1)./(abs(xest(:,tt)).^2 + ss + b);
    
            gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muzP(tt) +d);
            gammaOmegaD(tt) = gammaOmegaP(tt); 
    
            %% update st
%             st(2:end,tt) = b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt)./...
%                 (b^a/gamma(a).*exp((a-1).*lnrho(2:end,tt) - b.*rho(2:end,tt)).*lambdaIn(:,tt) + ...
%                 bbar^abar/gamma(abar).*exp((abar-1).*lnrho(2:end,tt) - bbar.*rho(2:end,tt)).*(1-lambdaIn(:,tt)));
            %% update D
            if (niterA == NiterA) && JESD 
                vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)' + Sx(:,:,tt))*Phid(:,:,tt)'));
                for ii = 1:dataNum
                    v(ii,tt)  = sum(vv(ii:dataNum:end));
                    mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
                    qq = -v(ii,tt)*gammaOmegaD(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2;
                    qq = qq - max(qq);
                    dPost(ii,tt) = sum(constell.*exp(qq))./sum(exp(qq));
                    sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(qq))./...
                        sum(exp(qq)) - abs(dPost(ii,tt))^2);
                end
    %                 dPost = dd;
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
%                 PhidPostSquared2(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
            end
        end
    end
end
%% VBI Iteration
for niter = 1:Niter
    xoldOut = xest;
    %% Subgraph A
    for niterA = 1:NiterA
        xold = xest;
        for tt = 1:T
            %% update D
            if (niter > 1 && niterA == 1) && JESD
                vv = abs(diag(Phid(:,:,tt)*(xest(:,tt)*xest(:,tt)' + Sx(:,:,tt))*Phid(:,:,tt)'));
                for ii = 1:dataNum
                    v(ii,tt)  = sum(vv(ii:dataNum:end));
                    mu(ii,tt) = conj(yd(ii:dataNum:end,tt)'*Phid(ii:dataNum:end,:,tt)*xest(:,tt));
        
                    qq = -v(ii,tt)*gammaOmegaD(tt)*abs(constell-mu(ii,tt)./v(ii,tt)).^2;
                    qq = qq - max(qq);
                    dPost(ii,tt) = sum(constell.*exp(qq))./sum(exp(qq));
                    sigmaDPost(ii,tt) = abs(sum(abs(constell).^2.*exp(qq))./...
                        sum(exp(qq)) - abs(dPost(ii,tt))^2);
                end
        %         PhidPost(:,:,tt) = sensingMatrixGenWSymbol(xUAV(:,tt),xR,dataSubcarrInd ,dPost(:,tt),0);
                PhidPost(:,:,tt) = repmat(kron(ones(NaNv,1),dPost(:,tt)),1,N).*Phid(:,:,tt);
        %         PhidPostSquared(:,:,tt) = PhidPost(:,:,tt)'*PhidPost(:,:,tt);
                PhidPostSquared(:,:,tt) = Phid(:,:,tt)'*(repmat(kron(ones(NaNv,1),abs(dPost(:,tt)).^2 + sigmaDPost(:,tt)),1,N).*Phid(:,:,tt));
            end
            %% update X
            Sx(:,:,tt) = inv((Phip(:,:,tt)'*Phip(:,:,tt))*gammaOmegaP(tt) +...
                PhidPostSquared(:,:,tt)*gammaOmegaD(tt) + diag(rho(:,tt)));
            xest(:,tt) = Sx(:,:,tt)*( Phip(:,:,tt)'*yp(:,tt)*gammaOmegaP(tt) +...
                PhidPost(:,:,tt)'*yd(:,tt)*gammaOmegaD(tt) );
            if (niter > 1) && JESD
                muz(tt) = abs((trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)')+...
                    xest(:,tt)'*(PhidPostSquared(:,:,tt)-PhidPost(:,:,tt)'*PhidPost(:,:,tt))*xest(:,tt)+...
                    trace(Sx(:,:,tt)*(PhidPostSquared(:,:,tt)) )));
                muzP(tt) = muz(tt);muzD(tt) = muz(tt);
            else
                muzP(tt) = abs(trace(Phip(:,:,tt)*Sx(:,:,tt)*Phip(:,:,tt)'));
            end
            ss = abs(diag(Sx(:,:,tt)));
            %% update rho and gammaOmega
            rho(:,tt) = (a + 1)./(abs(xest(:,tt)).^2 + ss + b);
            
            if (niter > 1) && JESD
%                 gammaOmegaP(tt) = (c+1)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2/Mp + muzP(tt) +d);
%                 gammaOmegaD(tt) = (c+1)./((norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2)/Md + muzD(tt) +d);

                gammaOmega(tt) = (c+(Md+Mp))./((norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2+...
                    norm(yd(:,tt) - PhidPost(:,:,tt)*xest(:,tt),2)^2) + muz(tt) +d);
                gammaOmegaP(tt) = gammaOmega(tt);gammaOmegaD(tt) = gammaOmega(tt);
            else
                gammaOmegaP(tt) = (c+Mp)./(norm(yp(:,tt) - Phip(:,:,tt)*xest(:,tt),2)^2 + muzP(tt) +d);
                gammaOmegaD(tt) = gammaOmegaP(tt); gammaOmega(tt) = gammaOmegaP(tt);
            end
            
        end
        %% Break when converging
        % to be completed
        nmse(niterA,niter) = mean(vecnorm(xest-gg,2).^2./vecnorm(gg,2).^2);
        if norm(xest - xold,2)^2/norm(xold)^2 < 1e-3
            break;
        end
        %
        %%%%%%%%%%%%%%%%%
    end
    %% EM update
    if JESD && ~convBreak
        dDecoded = mleDecoder(dPost,constell);
        SER(niter) = sum(dDecoded ~= dd,'all')/numel(dd);
    end
%     lambda0 = mean(st(2:end,:),'all');
    if norm(xest - xoldOut,'fro')^2/norm(xoldOut,'fro')^2 < 1e-3 && convBreak
        break;
    end
end
state.NMSE = nmse;
state.rho = rho;
state.gammaOmega = gammaOmega;
state.lambda0 = lambda0;
if JESD && ~convBreak
    state.SER = SER;
end
% state.telap = toc;
end