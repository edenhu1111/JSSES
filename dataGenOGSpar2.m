%% ISAC Tx & Rx
clear all;
Nsize = 20;numBS = 2;

fprintf('Data Generator for Sparsity\n');
load('dataset\3\PhipGenie.mat');
PhipGenie = PhipGeniePO{2,3};
sparArray = [1,3,5,7,9];
ind{1} = [5];
ind{2} = [2,4,5];
ind{3} = [1:2,4:5,9];
ind{4} = [1:2,4:8];
ind{5} = [1:9];
x = cell(length(sparArray),1);
zp = cell(length(sparArray),1);
xx = (5 + randn(9,numBS)).*exp(1j*2*pi*rand(9,numBS));
for ii = 1:length(sparArray)
    x{ii} = zeros(Nsize^2,numBS);
    indTmp = [(Nsize)/4+1,2*(Nsize)/4+1,3*(Nsize)/4+1];
    indTmp = round(indTmp);
    targetInd = [indTmp + Nsize*(round((Nsize)/4)),indTmp + 2*Nsize*(round((Nsize)/4)),indTmp + 3*Nsize*(round((Nsize)/4))];
    targetInd = targetInd(ind{ii});
    numTarget = length(targetInd);
    x{ii}(targetInd,:) =  xx(ind{ii},:);
    for jj = 1:numBS
        zp{ii}(:,jj) = PhipGenie(:,:,jj)*x{ii}(:,jj);
    end
end
save('dataset\1\sparArray','sparArray');
save('dataset\1\zpSpar2','zp');
save('dataset\1\xSpar2','x');

