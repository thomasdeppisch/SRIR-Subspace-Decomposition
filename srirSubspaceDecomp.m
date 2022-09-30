function [dirSrir, resSrir, numDirSubspaceComponents, gsvs, detectionThreshold, gsvSum, avgGsvSum] = ...
            srirSubspaceDecomp(srir, fs, blockLenSmp, hopSizeSmp, kappa, numBlocksGsvSumAvg, residualEstimateLengthMs, ...
                               decompositionTimeLimitMs, numBlocksSmoothThresh)
% Direct and Residual Subspace Decomposition of Spatial Room Impulse
% Responses, for details please refer to the following publication:
%
% T. Deppisch, S. Amengual Gari, P. Calamia, J. Ahrens, "Direct and
% Residual Subspace Decomposition of Spatial Room Impulse Responses,"
% arXiv preprint, 2022, doi:10.48550/arxiv.2207.09733.
%
% thomas.deppisch@chalmers.se

% Input parameters (many have default values)
arguments
    srir (:,:)                               % Spatial room impulse response, samples x channels
    fs (1,1)                                 % Sampling frequency in Hz
    blockLenSmp (1,1) = max(32,size(srir,2)) % Block length in samples
    hopSizeSmp (1,1) = round(blockLenSmp/8)  % Hop size in samples
    kappa (1,1) = 3                          % Detection threshold: kappa is the number of standard deviations that the sum of the generalized singular values need to exceed over their mean so that a reflection is detected.
    numBlocksGsvSumAvg (1,1) = 32            % Number of blocks to average the sum of the generalized singular values to determine the detection threshold.
    residualEstimateLengthMs (1,1) = 20      % Length of the residual (noise) estimate in milliseconds
    decompositionTimeLimitMs (1,1) = 100     % Latest time instant at which the decomposition is performed in milliseconds. (Time limit after which no salient reflections are expected.)
    numBlocksSmoothThresh (1,1) = 1          % Number of signal blocks over which the estimation of the number of subspace components is smoothed. (No smoothing if numBlocksSmoothThresh = 1. If numBlocksSmoothThresh = 3, the number of subspace components of the current block is the maximum of the number of estimated subspace components of the current and the two next blocks.)
end

% Output:
% dirSrir                  .. The direct part of the SRIR.
% resSrir                  .. The residual of the SRIR.
% numDirSubspaceComponents .. The number of direct subspace components used in the decomposition of each block.
% gsvs                     .. The generalized singular values of each block.
% detectionThreshold       .. The detection threshold. The subspace decomposition if performed if the sum of the GSVs exceeds the detection threshold.
% gsvSum                   .. The sum of the generalized singular values.
% avgGsvSum                .. The time-averaged sum of the GSVs that is used in the calculation of the detectionThreshold.

numSamples = size(srir,1);
numChannels = size(srir,2);

numValidSamples = hopSizeSmp;
numResidualSamples = residualEstimateLengthMs / 1000 * fs;
decompositionTimeLimitSmp = decompositionTimeLimitMs / 1000 * fs;

% some sanity checks
if mod(blockLenSmp,2) ~= 0 || mod(hopSizeSmp,2) ~= 0
    error('srirSubspaceDecomp: winLenSmp and stepSizeSmp must both be even integers!')
end
if blockLenSmp < numChannels
    error('srirSubspaceDecomp: winLenSmp needs to be equal or larger than numChannels.')
end
if numResidualSamples < numChannels
    error('srirSubspaceDecomp: numResidualSamples needs to be equal or larger than numChannels.')
end
if decompositionTimeLimitSmp + numResidualSamples > numSamples
    error('srirSubspaceDecomp: decompositionTimeLimitSmp plus numResidualSamples exceeds the impulse response.')
end

win = hann(blockLenSmp);
overlapSmp = blockLenSmp - hopSizeSmp;

srirPadded = [zeros(overlapSmp/2, numChannels); srir];
numBlocks = ceil((size(srirPadded,1) - overlapSmp) / hopSizeSmp);
numNoiseBlocks = round((numResidualSamples - blockLenSmp) / numValidSamples + 1);
srirBuffered = zeros(blockLenSmp, numChannels, numBlocks);
for ch = 1:numChannels
    srirBuffered(:,ch,:) = buffer(srirPadded(:,ch), blockLenSmp, overlapSmp, 'nodelay');
end

blockOffset = overlapSmp/2 + 1;
noiseDataMatrix = zeros(numResidualSamples, numChannels);
dirSrir = zeros(numSamples, numChannels);
resSrir = srir;
numComponents = zeros(numBlocks,1);
numSingVals = min(numChannels,blockLenSmp);
gsvs = zeros(numBlocks, numSingVals); % generalized singular values
singVals = zeros(numSingVals, numBlocks);
lSingVecs = zeros(blockLenSmp, numSingVals, numBlocks);
rSingVecs = zeros(numChannels, numChannels, numBlocks);
gsvSum = zeros(numBlocks, 1);
numDirSubspaceComponents = zeros(numBlocks, 1);
gsvSumBuffer = zeros(numBlocksGsvSumAvg,1);
avgGsvSum = zeros(numBlocks,1);
gsvStDevDuringNoise = zeros(numBlocks,1);
singValIdxOffset = max(0,numChannels-blockLenSmp);

for ii = numBlocks-numBlocksSmoothThresh+1:-1:1
    % idx at which subspace decomposition is performed (decomposition is
    % delayed if numBlocksSmoothThresh is larger than 1):
    decompBlockIdx = ii + numBlocksSmoothThresh - 1; % skewed window

    blockStartSmpIdx = (decompBlockIdx-1) * hopSizeSmp + 1;
    blockEndSmpIdx = blockStartSmpIdx + numValidSamples - 1;

    XWin = srirBuffered(:,:,ii) .* win;

    if ii <= numBlocks-numNoiseBlocks
        % perform the GSVD
        [lSingVecs(:,:,ii), ~, rSingVecs(:,:,ii), singValsTemp, noiseSingValsTemp] = gsvd(XWin, noiseDataMatrix, 0);
        singVals(:,ii) = diag(singValsTemp, singValIdxOffset);
        gsvs(ii,:) = singVals(:,ii).^2./diag(noiseSingValsTemp(singValIdxOffset+1:end, singValIdxOffset+1:end)).^2; % ordered increasingly
        gsvSum(ii) = sum(gsvs(ii,:));
    
        if blockStartSmpIdx <= decompositionTimeLimitSmp % do the actual subspace decomposition
            if gsvSum(ii) > (avgGsvSum(ii+1) + kappa * gsvStDevDuringNoise(ii+1)) % detection threshold
                numSingValsNoiseSubspace = sum(cumsum(gsvs(ii, :)) .* numSingVals./(1:numSingVals) < avgGsvSum(ii+1)); % number of direct subspace components threshold
                numComponents(ii) = numSingVals - numSingValsNoiseSubspace; % number of direct subspace components
            else
                numComponents(ii) = 0;
            end

            % perform a smoothing by taking the maximum of the numSources over a number of analysis blocks (this does nothing if numBlocksSmoothThresh = 1)
            numDirSubspaceComponents(decompBlockIdx) = max(numComponents(ii:ii+numBlocksSmoothThresh-1));
        
            % subspace decomposition
            U = lSingVecs(blockOffset:blockOffset+numValidSamples-1, :, decompBlockIdx);
        
            Ss = zeros(numSingVals, numChannels);
            Sn = zeros(numSingVals, numChannels);
            Ss(:, singValIdxOffset+1:end) = diag([zeros(numSingVals-numDirSubspaceComponents(decompBlockIdx),1); singVals(end-numDirSubspaceComponents(decompBlockIdx)+1:end, decompBlockIdx)]);
            Sn(:, singValIdxOffset+1:end) = diag([singVals(1:1:end-numDirSubspaceComponents(decompBlockIdx), decompBlockIdx); zeros(numDirSubspaceComponents(decompBlockIdx),1)]);
        
            % calculate direct and residual subspace signal blocks
            dirSrir(blockStartSmpIdx:blockEndSmpIdx, :) = U * Ss * rSingVecs(:,:,decompBlockIdx)';
            resSrir(blockStartSmpIdx:blockEndSmpIdx, :) = U * Sn * rSingVecs(:,:,decompBlockIdx)';
        end
    end

    % accumulate noise data if no subspace decomposition is performed
    if numDirSubspaceComponents(decompBlockIdx) == 0
        noiseDataMatrix = circshift(noiseDataMatrix, numValidSamples, 1);
        noiseDataMatrix(1:numValidSamples,:) = srirBuffered(blockOffset:blockOffset+numValidSamples-1, :, decompBlockIdx);
    end

    % save the sum of the GSVs in the noisy case (no salient reflection occuring)
    if ii <= numBlocks-numNoiseBlocks
        if numComponents(ii) == 0 % only update the sum if there are no reflections detected
            gsvSumBuffer = circshift(gsvSumBuffer,1);
            gsvSumBuffer(1) = sum(gsvs(ii,:));
            avgGsvSum(ii) = mean(gsvSumBuffer);
            gsvStDevDuringNoise(ii) = std(gsvSumBuffer);
        else
            avgGsvSum(ii) = avgGsvSum(ii+1);
            gsvStDevDuringNoise(ii) = gsvStDevDuringNoise(ii+1);
        end
    end

end

% output detection threshold for analysis purposes
detectionThreshold = avgGsvSum + kappa * gsvStDevDuringNoise;
