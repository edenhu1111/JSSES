function [outTheta] = CoarseEstimation(y,inputParam)
sMGParam = inputParam.sMGParam;
Na = sMGParam.Na;
y = (reshape(y,[],Na));

[S,wo] = pmusic(y,10);
S = fftshift(S)./max(S);
ind = find(S>2e-1);
[~,locs] = findpeaks(S);
locs = intersect(ind,locs);
outTheta = (wo - pi)/pi;
outTheta = outTheta(locs);
outTheta = acos(outTheta.');
end