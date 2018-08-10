% generate multiple networks with different rewiring
% probabilities and assess their small-world measures

% random seed
randomseed = randi(100);
rng(randomseed);

% number of nodes in networks
N = 250;
% 2 * number of connected neighbors (K on either side)
K = 25;
% rewiring probabilities
q = 0:0.05:0.2;
% number of networks to simulate in each condition ("samples")
numTrialNetworks = 100;
% keep track of statistics for these networks
smallWorldMeasure = NaN(length(q), numTrialNetworks);
clusterCoeffRatio = NaN(length(q), numTrialNetworks);   % gamma
charPathLengthRatio = NaN(length(q), numTrialNetworks); % lambda

for qix = 1:length(q)
    tic
    disp(['*** q = ' num2str(q(qix))]);
    for t = 1:numTrialNetworks
        if mod(t,10)==0
            disp(t);
        end
        networkObject = createNetwork(N, K, q(qix), 0, 0);
        smallWorldMeasure(qix,t) = networkObject.smallWorldMeasure;
        clusterCoeffRatio(qix,t) = networkObject.clusteringCoefficientRatio;
        charPathLengthRatio(qix,t) = networkObject.characteristicPathLengthRatio;
    end
    toc
end

% compute mean of all measures as a function of q
meanSW = mean(smallWorldMeasure,2);
meanCC = mean(clusterCoeffRatio,2);
meanPL = mean(charPathLengthRatio,2);

% compute error bars (standard error of the mean) for each rewiring probability
semSW = std(smallWorldMeasure,[],2) ./ sqrt(numTrialNetworks);
semCC = std(clusterCoeffRatio,[],2) ./ sqrt(numTrialNetworks);
semPL = std(smallWorldMeasure,[],2) ./ sqrt(numTrialNetworks);

% visualize
% small-world-measure
figure(1);
clf;
plot(q, meanSW);
hold on;
errorbar(q, meanSW, semSW, 'xr');
xlabel('Rewiring probability');
ylabel('Small-world measure');
title(['N = ' num2str(N) ', K = ' num2str(K) ', nsamples = ' num2str(numTrialNetworks)]);
hold off;
saveas(gca, ['smallWorldMeasure_N' num2str(N) '_K' num2str(K) '_' num2str(numTrialNetworks)]);

% cluster coefficient ratio (gamma) and characteristic path length ratio
% (lambda)
figure(2);
clf;
hold on;
plot(q, meanCC, 'displayname', 'cluster coefficient');
plot(q, meanPL, 'displayname', 'avg. path length');
legend('show');
errorbar(q, meanCC, semCC, 'xr');
errorbar(q, meanPL, semPL, 'xr');
xlabel('Rewiring probability');
title(['N = ' num2str(N) ', K = ' num2str(K) ', nsamples = ' num2str(numTrialNetworks)]);
hold off;
saveas(gca, ['gamma&lambda_N' num2str(N) '_K' num2str(K) '_' num2str(numTrialNetworks)]);


