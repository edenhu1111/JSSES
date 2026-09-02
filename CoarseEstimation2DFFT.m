function [angleDelayEst, state] = CoarseEstimation2DFFT(y, inputParam)
%COARSEESTIMATION2DFFT Estimate direction cosine and delay using a 2-D FFT.
%   ANGLEDELAYEST = COARSEESTIMATION2DFFT(Y, INPUTPARAM) reshapes the
%   antenna/subcarrier samples in Y, applies an oversampled two-dimensional
%   inverse FFT, and extracts separated peaks. Each row of ANGLEDELAYEST is
%       [directionCosine, totalDelaySeconds].
%
%   Required INPUTPARAM fields:
%       numAntennas        Number of ULA elements.
%       subcarrierIndices  Uniformly spaced OFDM subcarrier indices.
%       deltaF             Subcarrier spacing in Hz.
%
%   Optional INPUTPARAM fields:
%       numTargets          Number of peaks to extract. If omitted, an MDL
%                           estimate is used.
%       oversampling        Scalar or [delay, angle] factors (default 8).
%       directionCosineRange Search interval, default [-1, 1].
%       delayRange          Delay interval in seconds, default one
%                           unambiguous delay period.
%       pilotSymbols        Known pilots to remove before the FFT.
%       suppressionRadius   Radius [delay, angle] in FFT bins.
%       maxNumTargets       MDL search upper bound (default 16).
%       refinePeaks         Use local quadratic sub-bin interpolation
%                           (default true).
%
%   Y follows the sensing-matrix ordering used in this repository:
%   subcarrier samples vary fastest, followed by antenna elements.

if ~isnumeric(y) || isempty(y)
    error('CoarseEstimation2DFFT:InvalidObservation', ...
        'y must be a nonempty numeric array.');
end
if ~isstruct(inputParam)
    error('CoarseEstimation2DFFT:InvalidConfiguration', ...
        'inputParam must be a configuration structure.');
end
requiredFields = {'numAntennas', 'subcarrierIndices', 'deltaF'};
for fieldIndex = 1:numel(requiredFields)
    if ~isfield(inputParam, requiredFields{fieldIndex})
        error('CoarseEstimation2DFFT:MissingConfiguration', ...
            'inputParam.%s is required.', requiredFields{fieldIndex});
    end
end
inputParam = setDefault(inputParam, 'numTargets', NaN);
inputParam = setDefault(inputParam, 'oversampling', 8);
inputParam = setDefault(inputParam, 'directionCosineRange', [-1, 1]);
inputParam = setDefault(inputParam, 'delayRange', [0, Inf]);
inputParam = setDefault(inputParam, 'pilotSymbols', []);
inputParam = setDefault(inputParam, 'suppressionRadius', []);
inputParam = setDefault(inputParam, 'maxNumTargets', 16);
inputParam = setDefault(inputParam, 'refinePeaks', true);

validateattributes(inputParam.numAntennas, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'numAntennas');
validateattributes(inputParam.subcarrierIndices, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'}, mfilename, 'subcarrierIndices');
validateattributes(inputParam.deltaF, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'deltaF');
validateattributes(inputParam.oversampling, {'numeric'}, ...
    {'vector', 'real', 'finite', 'positive'}, mfilename, 'oversampling');
validateattributes(inputParam.directionCosineRange, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite'}, mfilename, ...
    'directionCosineRange');
validateattributes(inputParam.delayRange, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'nonnegative'}, mfilename, 'delayRange');

numAntennas = inputParam.numAntennas;
subcarrierIndices = inputParam.subcarrierIndices(:);
numSubcarriers = numel(subcarrierIndices);

if numel(y) ~= numSubcarriers * numAntennas
    error('CoarseEstimation2DFFT:InvalidInputSize', ...
        ['Expected %d samples (%d subcarriers x %d antennas), but ', ...
         'received %d.'], numSubcarriers * numAntennas, ...
        numSubcarriers, numAntennas, numel(y));
end

subcarrierSteps = diff(subcarrierIndices);
if isempty(subcarrierSteps)
    subcarrierStride = 1;
elseif any(abs(subcarrierSteps - subcarrierSteps(1)) > ...
        10 * eps(max(abs(subcarrierIndices))))
    error('CoarseEstimation2DFFT:NonuniformSubcarriers', ...
        'The 2-D FFT estimator requires uniformly spaced subcarriers.');
else
    subcarrierStride = subcarrierSteps(1);
end
if subcarrierStride <= 0
    error('CoarseEstimation2DFFT:InvalidSubcarrierOrder', ...
        'subcarrierIndices must be strictly increasing.');
end

sampleMatrix = reshape(y(:), numSubcarriers, numAntennas);
if ~isempty(inputParam.pilotSymbols)
    pilots = inputParam.pilotSymbols;
    if isvector(pilots) && numel(pilots) == numSubcarriers
        pilots = repmat(pilots(:), 1, numAntennas);
    elseif ~isequal(size(pilots), size(sampleMatrix))
        error('CoarseEstimation2DFFT:InvalidPilotSize', ...
            'pilotSymbols must have one entry per subcarrier or per sample.');
    end
    if any(abs(pilots(:)) < sqrt(eps))
        error('CoarseEstimation2DFFT:ZeroPilot', ...
            'pilotSymbols cannot contain zeros.');
    end
    sampleMatrix = sampleMatrix ./ pilots;
end

oversampling = inputParam.oversampling;
if isscalar(oversampling)
    oversampling = [oversampling, oversampling];
elseif numel(oversampling) ~= 2
    error('CoarseEstimation2DFFT:InvalidOversampling', ...
        'oversampling must be a scalar or [delay, angle].');
end
fftSize = max([numSubcarriers, numAntennas], ...
    ceil([numSubcarriers, numAntennas] .* oversampling));

% Both model dimensions contain negative-exponent complex sinusoids, hence
% IFFT is used along both dimensions. Only the signed angle dimension is
% shifted; physical delay is restricted to its nonnegative unambiguous bin.
spectrumComplex = ifft(sampleMatrix, fftSize(1), 1);
spectrumComplex = ifft(spectrumComplex, fftSize(2), 2);
spectrumComplex = fftshift(spectrumComplex, 2);
spectrumPower = abs(spectrumComplex).^2;

delayFrequency = (0:fftSize(1)-1).' / fftSize(1);
directionCosine = 2 * ((0:fftSize(2)-1) - floor(fftSize(2)/2)) / fftSize(2);
delayAxis = delayFrequency / (subcarrierStride * inputParam.deltaF);

directionRange = sort(inputParam.directionCosineRange);
delayRange = sort(inputParam.delayRange);
if isinf(delayRange(2))
    delayRange(2) = delayAxis(end);
end
validDelay = delayAxis >= max(0, delayRange(1)) & delayAxis <= delayRange(2);
validAngle = directionCosine >= max(-1, directionRange(1)) & ...
    directionCosine <= min(1, directionRange(2));
searchMask = validDelay * validAngle;
spectrumForSearch = spectrumPower;
spectrumForSearch(~logical(searchMask)) = -Inf;

if isnan(inputParam.numTargets)
    numTargets = estimateTargetCount(sampleMatrix, inputParam.maxNumTargets);
else
    numTargets = inputParam.numTargets;
end
numTargets = min(numTargets, nnz(searchMask));

if isempty(inputParam.suppressionRadius)
    suppressionRadius = max(1, ceil(fftSize ./ ...
        [numSubcarriers, numAntennas]));
else
    suppressionRadius = inputParam.suppressionRadius;
    if isscalar(suppressionRadius)
        suppressionRadius = [suppressionRadius, suppressionRadius];
    elseif numel(suppressionRadius) ~= 2
        error('CoarseEstimation2DFFT:InvalidSuppressionRadius', ...
            'suppressionRadius must be a scalar or [delay, angle].');
    end
end
suppressionRadius = ceil(suppressionRadius);

peakBins = zeros(numTargets, 2);
peakPower = zeros(numTargets, 1);
workSpectrum = spectrumForSearch;
numFound = 0;
for targetIdx = 1:numTargets
    [currentPower, linearIndex] = max(workSpectrum(:));
    if ~isfinite(currentPower)
        break;
    end
    numFound = numFound + 1;
    [delayBin, angleBin] = ind2sub(fftSize, linearIndex);
    peakBins(numFound, :) = [delayBin, angleBin];
    peakPower(numFound) = currentPower;

    delayBins = max(1, delayBin-suppressionRadius(1)): ...
        min(fftSize(1), delayBin+suppressionRadius(1));
    angleBins = mod((angleBin-suppressionRadius(2): ...
        angleBin+suppressionRadius(2))-1, fftSize(2)) + 1;
    workSpectrum(delayBins, angleBins) = -Inf;
end

peakBins = peakBins(1:numFound, :);
peakPower = peakPower(1:numFound);
peakOffsets = zeros(numFound, 2);
angleDelayEst = zeros(numFound, 2);
if numFound > 0
    if inputParam.refinePeaks
        logSpectrum = log(max(spectrumPower, realmin));
        for peakIdx = 1:numFound
            peakOffsets(peakIdx, :) = quadraticPeakOffset( ...
                logSpectrum, peakBins(peakIdx, :));
        end
    end
    continuousDelayBin = peakBins(:,1)-1 + peakOffsets(:,1);
    continuousAngleBin = peakBins(:,2)-1 + peakOffsets(:,2);
    angleDelayEst(:,1) = 2 * (continuousAngleBin - ...
        floor(fftSize(2)/2)) / fftSize(2);
    angleDelayEst(:,2) = continuousDelayBin / fftSize(1) / ...
        (subcarrierStride * inputParam.deltaF);
    [~, order] = sort(peakPower, 'descend');
    angleDelayEst = angleDelayEst(order, :);
    peakBins = peakBins(order, :);
    peakOffsets = peakOffsets(order, :);
    peakPower = peakPower(order);
end

state.fftSize = fftSize;
state.directionCosineAxis = directionCosine;
state.delayAxis = delayAxis;
state.spectrumPower = spectrumPower;
state.peakBins = peakBins;
state.peakOffsets = peakOffsets;
state.peakPower = peakPower;
state.numTargets = numFound;
state.unambiguousDelay = 1 / (subcarrierStride * inputParam.deltaF);
end

function offsets = quadraticPeakOffset(logSpectrum, peakBin)
% Independent one-dimensional parabolic fits through the local 2-D peak.
numDelayBins = size(logSpectrum, 1);
numAngleBins = size(logSpectrum, 2);
delayBin = peakBin(1);
angleBin = peakBin(2);
offsets = [0, 0];

if delayBin > 1 && delayBin < numDelayBins
    values = logSpectrum(delayBin-1:delayBin+1, angleBin);
    denominator = values(1) - 2*values(2) + values(3);
    if isfinite(denominator) && abs(denominator) > eps
        offsets(1) = 0.5 * (values(1)-values(3)) / denominator;
    end
end

previousAngle = mod(angleBin-2, numAngleBins) + 1;
nextAngle = mod(angleBin, numAngleBins) + 1;
values = logSpectrum(delayBin, [previousAngle, angleBin, nextAngle]);
denominator = values(1) - 2*values(2) + values(3);
if isfinite(denominator) && abs(denominator) > eps
    offsets(2) = 0.5 * (values(1)-values(3)) / denominator;
end
offsets = max(-0.5, min(0.5, offsets));
end

function numTargets = estimateTargetCount(sampleMatrix, maxNumTargets)
% Wax-Kailath MDL estimate applied to the antenna/subcarrier data matrix.
[numRows, numColumns] = size(sampleMatrix);
numSensors = min(numRows, numColumns);
numSnapshots = max(numRows, numColumns);
singularValues = svd(sampleMatrix, 'econ');
eigenvalues = singularValues(1:numSensors).^2 / numSnapshots;
maxOrder = min([maxNumTargets, numSensors-2]);
mdl = inf(maxOrder+1, 1);
for order = 0:maxOrder
    noiseEigenvalues = max(eigenvalues(order+1:end), realmin);
    geometricMean = exp(mean(log(noiseEigenvalues)));
    arithmeticMean = mean(noiseEigenvalues);
    likelihoodTerm = -numSnapshots * (numSensors-order) * ...
        log(max(geometricMean/arithmeticMean, realmin));
    penaltyTerm = 0.5 * order * (2*numSensors-order) * log(numSnapshots);
    mdl(order+1) = likelihoodTerm + penaltyTerm;
end
[~, minimumIndex] = min(mdl);
numTargets = minimumIndex - 1;
end

function outputStruct = setDefault(inputStruct, fieldName, defaultValue)
outputStruct = inputStruct;
if ~isfield(outputStruct, fieldName) || isempty(outputStruct.(fieldName))
    outputStruct.(fieldName) = defaultValue;
end
end
