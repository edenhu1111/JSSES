function [pTxRefined, biasRangeRefined, finalTargets, state] = ...
    estimate_bias_by_coherence(pTxPrior, pRx, theta1, range1, ...
    theta2, range2, options)
%ESTIMATE_BIAS_BY_COHERENCE Associate two-Rx peaks and solve by adaptive LM.
%   RANGE1 and RANGE2 are total bistatic path lengths in metres, including
%   the common clock-bias range c*DeltaTau. THETA1 and THETA2 are global
%   azimuth angles in radians.
%
%   The function first intersects AoA rays and associates the two-Rx
%   measurements using their clock- and transmitter-independent range
%   difference. It then solves for [pTxX, pTxY, c*DeltaTau] using a
%   Levenberg-Marquardt iteration with an adaptive damping parameter.

if nargin < 7 || isempty(options)
    options = struct();
end
options = setDefault(options, 'associationThreshold', 30);
options = setDefault(options, 'parallelTolerance', 0.05);
options = setDefault(options, 'maxIterations', 50);
options = setDefault(options, 'maxDampingTrials', 10);
options = setDefault(options, 'stepTolerance', 1e-6);
options = setDefault(options, 'costTolerance', 1e-10);
options = setDefault(options, 'initialLambda', 1e-2);
options = setDefault(options, 'lambdaIncrease', 10);
options = setDefault(options, 'lambdaDecrease', 3);

pTxPrior = reshape(pTxPrior, 1, []);
if numel(pTxPrior) ~= 2 || size(pRx, 2) ~= 2 || size(pRx, 1) ~= 2
    error('estimate_bias_by_coherence:InvalidGeometry', ...
        'pTxPrior must be 1x2 and pRx must be 2x2.');
end
theta1 = theta1(:);
theta2 = theta2(:);
range1 = range1(:);
range2 = range2(:);
if numel(theta1) ~= numel(range1) || numel(theta2) ~= numel(range2)
    error('estimate_bias_by_coherence:InvalidMeasurements', ...
        'Each angle vector must have the same length as its range vector.');
end

numPeaks1 = numel(theta1);
numPeaks2 = numel(theta2);
receiver1 = pRx(1, :);
receiver2 = pRx(2, :);
associationCost = inf(numPeaks1, numPeaks2);
candidatePositions = cell(numPeaks1, numPeaks2);

for peak1 = 1:numPeaks1
    direction1 = [cos(theta1(peak1)); sin(theta1(peak1))];
    for peak2 = 1:numPeaks2
        direction2 = [cos(theta2(peak2)); sin(theta2(peak2))];
        rayMatrix = [direction1, -direction2];
        if abs(det(rayMatrix)) <= options.parallelTolerance
            continue;
        end
        rayDistances = rayMatrix \ (receiver2 - receiver1).';
        if any(rayDistances <= 0)
            continue;
        end
        targetPosition = receiver1 + rayDistances(1) * direction1.';
        candidatePositions{peak1, peak2} = targetPosition;

        % The transmitter-to-target distance and common clock bias cancel
        % in this difference, making it the appropriate association cost.
        measuredDifference = range1(peak1) - range2(peak2);
        geometricDifference = rayDistances(1) - rayDistances(2);
        differenceResidual = abs(measuredDifference - geometricDifference);
        if differenceResidual <= options.associationThreshold
            associationCost(peak1, peak2) = differenceResidual;
        end
    end
end

if isempty(associationCost) || all(~isfinite(associationCost), 'all')
    error('estimate_bias_by_coherence:NoAssociation', ...
        'No physically valid two-Rx angle-delay association was found.');
end
matches = matchpairs(associationCost, options.associationThreshold);
if isempty(matches)
    error('estimate_bias_by_coherence:NoAssociation', ...
        'No two-Rx measurement pair passed the association threshold.');
end

numTargets = size(matches, 1);
finalTargets = zeros(numTargets, 2);
observedRanges = zeros(2*numTargets, 1);
for targetIdx = 1:numTargets
    peak1 = matches(targetIdx, 1);
    peak2 = matches(targetIdx, 2);
    finalTargets(targetIdx, :) = candidatePositions{peak1, peak2};
    observedRanges(2*targetIdx-1:2*targetIdx) = ...
        [range1(peak1); range2(peak2)];
end

state.matches = matches;
state.associationCost = associationCost;
initialParameter = [pTxPrior(:); 0];
[initialResidual, initialJacobian] = buildResidualAndJacobian( ...
    initialParameter, finalTargets, pRx, observedRanges);
state.jacobianRank = rank(initialJacobian);
state.identifiable = state.jacobianRank == 3;

% After target positions have been fixed by AoA intersections, both Rx
% measurements of one target contain the same transmitter-distance-plus-
% bias equation. At least three geometrically diverse targets are therefore
% needed for three independent equations in [pTxX, pTxY, biasRange].
if ~state.identifiable
    pTxRefined = pTxPrior;
    biasRangeRefined = mean(initialResidual);
    state.exitReason = 'rank_deficient_geometry';
    state.iterations = 0;
    state.costHistory = [];
    state.lambdaHistory = [];
    state.finalResidual = initialResidual - biasRangeRefined;
    return;
end

parameter = initialParameter;
lambda = options.initialLambda;
costHistory = nan(options.maxIterations+1, 1);
lambdaHistory = nan(options.maxIterations, 1);
residual = initialResidual;
jacobian = initialJacobian;
cost = 0.5 * (residual' * residual);
costHistory(1) = cost;
exitReason = 'maximum_iterations';

for iteration = 1:options.maxIterations
    normalMatrix = jacobian' * jacobian;
    diagonalScale = max(diag(normalMatrix), sqrt(eps));
    accepted = false;

    for dampingTrial = 1:options.maxDampingTrials
        step = (normalMatrix + lambda*diag(diagonalScale)) \ ...
            (jacobian' * residual);
        trialParameter = parameter + step;
        trialResidual = buildResidualAndJacobian(trialParameter, ...
            finalTargets, pRx, observedRanges);
        trialCost = 0.5 * (trialResidual' * trialResidual);

        if trialCost < cost
            previousCost = cost;
            parameter = trialParameter;
            residual = trialResidual;
            cost = trialCost;
            lambda = max(lambda/options.lambdaDecrease, eps);
            accepted = true;
            break;
        end
        lambda = min(lambda*options.lambdaIncrease, 1/eps);
    end

    lambdaHistory(iteration) = lambda;
    costHistory(iteration+1) = cost;
    if ~accepted
        exitReason = 'no_descent_step';
        break;
    end
    if norm(step) <= options.stepTolerance * ...
            (norm(parameter) + options.stepTolerance)
        exitReason = 'step_tolerance';
        break;
    end
    if previousCost-cost <= options.costTolerance * max(1, previousCost)
        exitReason = 'cost_tolerance';
        break;
    end
    [residual, jacobian] = buildResidualAndJacobian(parameter, ...
        finalTargets, pRx, observedRanges);
end

pTxRefined = parameter(1:2).';
biasRangeRefined = parameter(3);
state.exitReason = exitReason;
state.iterations = iteration;
state.costHistory = costHistory(1:iteration+1);
state.lambdaHistory = lambdaHistory(1:iteration);
state.finalResidual = residual;
end

function [residual, jacobian] = buildResidualAndJacobian(parameter, ...
    targets, receivers, observations)
numTargets = size(targets, 1);
residual = zeros(2*numTargets, 1);
jacobian = zeros(2*numTargets, 3);
transmitterPosition = parameter(1:2).';
biasRange = parameter(3);

for targetIdx = 1:numTargets
    transmitterVector = transmitterPosition - targets(targetIdx, :);
    transmitterDistance = max(norm(transmitterVector), sqrt(eps));
    transmitterGradient = transmitterVector / transmitterDistance;
    for receiverIdx = 1:2
        rowIndex = 2*targetIdx - 2 + receiverIdx;
        receiverDistance = norm(targets(targetIdx,:) - ...
            receivers(receiverIdx,:));
        predictedRange = transmitterDistance + receiverDistance + biasRange;
        residual(rowIndex) = observations(rowIndex) - predictedRange;
        jacobian(rowIndex, :) = [transmitterGradient, 1];
    end
end
end

function outputStruct = setDefault(inputStruct, fieldName, defaultValue)
outputStruct = inputStruct;
if ~isfield(outputStruct, fieldName) || isempty(outputStruct.(fieldName))
    outputStruct.(fieldName) = defaultValue;
end
end
