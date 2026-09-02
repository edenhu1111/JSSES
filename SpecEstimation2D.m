function [centroids, state] = SpecEstimation2D(y, inputPhase1)
%SPECESTIMATION2D Compatibility wrapper for the 2-D FFT coarse estimator.
%   The legacy output uses [u, 2*nu], where u is direction cosine and
%   nu = subcarrierStride*deltaF*delay. New code should call
%   CoarseEstimation2DFFT directly to obtain delay in seconds.

requiredFields = {'numAntennas', 'subcarrierIndices', 'deltaF'};
for fieldIndex = 1:numel(requiredFields)
    if ~isfield(inputPhase1, requiredFields{fieldIndex})
        error('SpecEstimation2D:MissingConfiguration', ...
            'inputPhase1.%s is required by the FFT estimator.', ...
            requiredFields{fieldIndex});
    end
end

fftInput = inputPhase1;
if isfield(inputPhase1, 'grid') && ~isempty(inputPhase1.grid)
    fftInput.directionCosineRange = ...
        [min(inputPhase1.grid(1,:)), max(inputPhase1.grid(1,:))];
    subcarrierStride = median(diff(inputPhase1.subcarrierIndices));
    fftInput.delayRange = ...
        [min(inputPhase1.grid(2,:)), max(inputPhase1.grid(2,:))] ./ ...
        (2 * subcarrierStride * inputPhase1.deltaF);
end

[angleDelayEst, state] = CoarseEstimation2DFFT(y, fftInput);
subcarrierStride = median(diff(inputPhase1.subcarrierIndices));
centroids = [angleDelayEst(:,1), ...
    2 * subcarrierStride * inputPhase1.deltaF * angleDelayEst(:,2)];
end
