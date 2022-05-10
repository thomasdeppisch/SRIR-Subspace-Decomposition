function [dirSrir, diffSrir, smoothNumSources, gsvs, detectionThreshold, gsvSum, avgGsvSumDuringNoise] = srirSubspaceDecomp(srir, fs, winLenSmp, ...
    stepSizeSmp, noiseAnalysisTimeMs, mixingTimeMs, numBlocksSmoothThresh, numBlocksSingValsAvg, evalThresholdMarginFactor)

numSamples = size(srir,1);
numChannels = size(srir,2);

numValidSamples = stepSizeSmp;
numNoiseSamples = noiseAnalysisTimeMs / 1000 * fs;
mixingTimeSmp = mixingTimeMs/1000 * fs;

if nargin < 8 || isempty(numBlocksSingValsAvg)
    numBlocksSingValsAvg = 32;
end
if nargin < 9 || isempty(evalThresholdMarginFactor)
    evalThresholdMarginFactor = 4;
end

if mod(winLenSmp,2) ~= 0 || mod(stepSizeSmp,2) ~= 0
    error('winLenSmp and stepSizeSmp must both be even integers!')
end

if mod(numBlocksSmoothThresh,2) == 0
    warning('numBlocksSmoothThresh should be odd for a symmetric analysis window')
end

if winLenSmp+numNoiseSamples < numChannels
    error('winLenSmp + numNoiseSamples needs to larger than numChannels')
end

if (mixingTimeSmp + numNoiseSamples) > numSamples
    error('mixing time plus noise data estimation length exceeds the impulse response')
end

win = hann(winLenSmp);
overlapSmp = winLenSmp - stepSizeSmp;

srirPadded = [zeros(overlapSmp/2, numChannels); srir];

numBlocks = ceil((size(srirPadded,1) - overlapSmp) / stepSizeSmp);
numNoiseBlocks = round((numNoiseSamples - winLenSmp) / numValidSamples + 1);
srirBuffered = zeros(winLenSmp, numChannels, numBlocks);
for ch = 1:numChannels
    srirBuffered(:,ch,:) = buffer(srirPadded(:,ch), winLenSmp, overlapSmp, 'nodelay');
end

blockOffset = overlapSmp/2 + 1;
noiseDataMatrix = zeros(numNoiseSamples, numChannels);
dirSrir = zeros(numSamples, numChannels);
diffSrir = srir;
numSources = zeros(numBlocks,1);

numSingVals = min(numChannels,winLenSmp);
gsvs = zeros(numBlocks, numSingVals); % generalized singular values
singVals = zeros(numSingVals, numBlocks);
lSingVecs = zeros(winLenSmp, numSingVals, numBlocks);
rSingVecs = zeros(numChannels, numChannels, numBlocks);
gsvSum = zeros(numBlocks, 1);

% the threshold will be regarded as exceeded if any sample within the
% smoothThreshWin exceeds the threshold
smoothNumSources = zeros(numBlocks, 1);

gsvSumBuffer = zeros(numBlocksSingValsAvg,1);
avgGsvSumDuringNoise = zeros(numBlocks,1);
gsvStDevDuringNoise = zeros(numBlocks,1);
singValIdxOffset = max(0,numChannels-winLenSmp);

for ii = numBlocks-numBlocksSmoothThresh+1:-1:1

    % idx at which subspace decomposition is performed (decomposition is
    % delayed if numBlocksSmoothThresh is larger than 1):
    % decompBlockIdx = ii + round(numBlocksSmoothThresh/2) - 1; % symmetric window
    decompBlockIdx = ii + numBlocksSmoothThresh - 1; % skewed window

    blockStartSmpIdx = (decompBlockIdx-1) * stepSizeSmp + 1;
    blockEndSmpIdx = blockStartSmpIdx + numValidSamples - 1;

    XWin = srirBuffered(:,:,ii) .* win;

    if ii <= numBlocks-numNoiseBlocks
        % do the GSVD to estimate a singular value threshold but only do
        % the subspace decomposition if ii < mixingTimeSmp

        % perform the GSVD
        [lSingVecs(:,:,ii), ~, rSingVecs(:,:,ii), singValsTemp, noiseSingValsTemp] = gsvd(XWin, noiseDataMatrix, 0);
        singVals(:,ii) = diag(singValsTemp, singValIdxOffset);
        gsvs(ii,:) = singVals(:,ii).^2./diag(noiseSingValsTemp(singValIdxOffset+1:end, singValIdxOffset+1:end)).^2; % ordered increasingly
        gsvSum(ii) = sum(gsvs(ii,:));
    
        if blockStartSmpIdx <= mixingTimeSmp % do the actual subspace decomposition
            if gsvSum(ii) > (avgGsvSumDuringNoise(ii+1) + evalThresholdMarginFactor * gsvStDevDuringNoise(ii+1)) % detection threshold
                numSingValsNoiseSubspace = sum(cumsum(gsvs(ii, :)) .* numSingVals./(1:numSingVals) < avgGsvSumDuringNoise(ii+1)); % number of direct subspace components threshold
                %numSingValsNoiseSubspace = sum(cumsum(gsvs(ii, :)) < avgGsvSumDuringNoise(ii+1));
                numSources(ii) = numSingVals - numSingValsNoiseSubspace; % number of direct subspace components
            else
                numSources(ii) = 0;
            end

            % perform a smoothing by taking the maximum of the numSources over a number of analysis blocks (this does nothing if numBlocksSmoothThresh = 1)
            smoothNumSources(decompBlockIdx) = max(numSources(ii:ii+numBlocksSmoothThresh-1));
        
            % subspace decomposition
            U = lSingVecs(blockOffset:blockOffset+numValidSamples-1, :, decompBlockIdx);
        
            Ss = zeros(numSingVals, numChannels);
            Sn = zeros(numSingVals, numChannels);
            Ss(:, singValIdxOffset+1:end) = diag([zeros(numSingVals-smoothNumSources(decompBlockIdx),1); singVals(end-smoothNumSources(decompBlockIdx)+1:end, decompBlockIdx)]);
            Sn(:, singValIdxOffset+1:end) = diag([singVals(1:1:end-smoothNumSources(decompBlockIdx), decompBlockIdx); zeros(smoothNumSources(decompBlockIdx),1)]);
        
            % calculate direct and diffuse subspace signal blocks
            dirSrir(blockStartSmpIdx:blockEndSmpIdx, :) = U * Ss * rSingVecs(:,:,decompBlockIdx)';
            diffSrir(blockStartSmpIdx:blockEndSmpIdx, :) = U * Sn * rSingVecs(:,:,decompBlockIdx)';
        end
    end

    % accumulate noise data if no subspace decomposition is done (if ii > numBlocks - numNoiseBlocks, we do nothing else)
    if smoothNumSources(decompBlockIdx) == 0
        noiseDataMatrix = circshift(noiseDataMatrix, numValidSamples, 1);
        noiseDataMatrix(1:numValidSamples,:) = srirBuffered(blockOffset:blockOffset+numValidSamples-1, :, decompBlockIdx);
    end

    % save the sum of the eigenvalues in the noisy case (no reflection occuring)
    if ii <= numBlocks-numNoiseBlocks
        if numSources(ii) == 0 % only update the sum if there are no reflections detected
            gsvSumBuffer = circshift(gsvSumBuffer,1);
            gsvSumBuffer(1) = sum(gsvs(ii,:));
            avgGsvSumDuringNoise(ii) = mean(gsvSumBuffer);
            gsvStDevDuringNoise(ii) = std(gsvSumBuffer);
        else
            avgGsvSumDuringNoise(ii) = avgGsvSumDuringNoise(ii+1);
            gsvStDevDuringNoise(ii) = gsvStDevDuringNoise(ii+1);
        end
    end

end

% output detection threshold for analysis purposes
detectionThreshold = avgGsvSumDuringNoise + evalThresholdMarginFactor * gsvStDevDuringNoise;
