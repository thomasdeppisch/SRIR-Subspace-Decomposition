clear all
close all

%% test the direct and residual subspace decomposition of an SRIR
srirStruct = load('resources/TA_lecture_hall_em32_ACN_N3D.mat');
srir = srirStruct.srir;
fs = srirStruct.fs;

% parameters for the subspace decomposition, see the function header of srirSubspaceDecomp for details
blockLenSmp = 32;
hopSizeSmp = blockLenSmp / 8;
kappa = 3;
numBlocksGsvSumAvg = 32;
residualEstimateLengthMs = 20;
decompositionTimeLimitMs = 100;
numBlocksSmoothThresh = 1;

[dirSrir, resSrir, numDirSubspaceComponents, gsvs, detectionThreshold, gsvSum, avgGsvSum] = ...
            srirSubspaceDecomp(srir, fs, blockLenSmp, hopSizeSmp, kappa, numBlocksGsvSumAvg, residualEstimateLengthMs, ...
                               decompositionTimeLimitMs, numBlocksSmoothThresh);


%% plot
srirLen = size(srir,1);
t = linspace(0, srirLen/fs-1/fs, srirLen).';
tBlocks = (0:hopSizeSmp:size(gsvs,1)*hopSizeSmp-hopSizeSmp)/fs;

figure
hold on
plot(t*1000, db(sum(abs(dirSrir), 2)), 'LineWidth', 2)
plot(t*1000, db(sum(abs(resSrir), 2)), 'LineWidth', 2)
xlim([0,100])
xlabel('$t$ (ms)')
ylabel('$\| \cdot \|$ (dB)')
legend('$\mathbf{x}_\mathrm{d}(t)$', '$\mathbf{x}_\mathrm{r}(t)$')
grid on

numChannels = size(srir,2);
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
legend([hGSVCumsum, hAvgGSVs, hDetectionThresh], 'cumulative sums of GSVs', 'subspace component threshold', 'detection threshold')
