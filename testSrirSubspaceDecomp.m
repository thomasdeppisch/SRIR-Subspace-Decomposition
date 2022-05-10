clear all
close all

%% test the direct and residual subspace decomposition of an SRIR
rirStruct = load('TA_lecture_hall_em32_Genelec_S1_R1.mat');
srir = rirStruct.irs;
fs = rirStruct.fs;

shOrder = 4;
% TODO: remove dependency on this function
srirShdFiltered = encodeMicArray(srir, fs, true, shOrder=shOrder);

% parameters for the subspace decomposition, see the function header for details:
% TODO: explain what the parameters do and how they behave
blockLenSmp = 32;
hopSizeSmp = blockLenSmp / 8;
noiseAnalysisTimeMs = 20;
decompositionTimeLimitMs = 100;
numBlocksSmoothThresh = 1;
numBlocksGsvSumAvg = 32;
kappa = 3;

[dirSrir, resSrir, numDirSubspaceComponents, gsvs, detectionThreshold, gsvSum, avgGsvSum] = ...
            srirSubspaceDecomp(srirShdFiltered, fs, blockLenSmp, hopSizeSmp, noiseAnalysisTimeMs, decompositionTimeLimitMs, ...
            numBlocksSmoothThresh, numBlocksGsvSumAvg, kappa);

%% plot
srirLen = size(srirShdFiltered,1);
t = linspace(0, srirLen/fs-1/fs, srirLen).';
tBlocks = (0:hopSizeSmp:size(gsvs,1)*hopSizeSmp-hopSizeSmp)/fs;

figure
hold on
plot(t*1000, db(sum(abs(dirSrir), 2)), 'LineWidth', 2)
plot(t*1000, db(sum(abs(resSrir), 2)), 'LineWidth', 2)
xlim([0,100])
%ylim(yLimitsDb)
xlabel('$t$ (ms)')
ylabel('$\| \cdot \|$ (dB)')
legend('$\mathbf{x}_\mathrm{d}(t)$', '$\mathbf{x}_\mathrm{r}(t)$')
grid on

numChannels = size(srirShdFiltered,2);
cumsumGSVsColors = copper(numChannels);
cumsumGSVs = cumsum(gsvs,2);

figure
hold on
for ii = 1:numChannels
    if ii == numChannels
        hGSVCumsum = plot(tBlocks*1000, cumsumGSVs(:,ii), 'Color', cumsumGSVsColors(ii,:), 'LineWidth', 1.5);
    else
        plot(tBlocks*1000, cumsumGSVs(:,ii)*numChannels/ii, 'Color', cumsumGSVsColors(ii,:), 'LineWidth', 1.5)
    end
end
hAvgGSVs = plot(tBlocks*1000, avgGsvSum, 'k:');
hDetectionThresh = plot(tBlocks*1000, detectionThreshold, 'k');
grid on
xlabel('$t$ (ms)')
xlim([0,100])
ylim([0,10])
legend([hGSVCumsum, hAvgGSVs, hDetectionThresh], 'cumulative sums of GSVs', 'detection threshold', 'subspace component threshold')
