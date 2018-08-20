%%
rep = 5;
maxstep = max(decision(rep, :, 2));
maxstep = maxstep;
if maxstep == 0
    maxstep = maxsteps;
end
figure;
subplot(1,2,1)
plot(1:maxstep, r(rep, :, 1:maxstep))
subplot(1,2,2)
plot(1:maxstep, lams(rep, :, 1:maxstep))

%% DT distributions (Habiba)
set(0, 'defaultAxesFontSize', 14);
nBins = 150;
rtNonZero = rt;
rtNonZero(rtNonZero == 0) = nan; %maxsteps+1;
figure('Color', 'w', 'Position', [800 100 1000 500]);
subplot(1,2,1)
[rtCounts, rtCenters] = hist(rtNonZero, nBins);
totalCount = sum(cumsum(rtCounts), 2);
xCutoff = find(totalCount > totalCount(end) - 1, 1);
bar(rtCenters, rtCounts / totalCount(end), 1)
xlim([0 xCutoff])
ylabel('p(t_{D} = t)')
xlabel('Timestep')
title('PDF');
winFracs = wins / sum(wins);
undecidedFrac = 1 - sum(wins) / nReps;
annotTxt = sprintf(' f_{L} = %0.3f\n f_{R} = %0.3f\n f_{?} = %0.3f', winFracs(1), winFracs(2), undecidedFrac);
annot = annotation('textbox', [.35 .6 .3 .3], 'String', annotTxt, 'FitBoxToText', 'on');
annot.FontSize = 12;
subplot(1,2,2)
cdf1 = cdfplot(rtNonZero(:, 1));
set(cdf1, 'LineWidth', 2);
hold on
cdf2 = cdfplot(rtNonZero(:, 2));
%xlim([0 maxsteps])
set(cdf2, 'LineWidth', 2);
set(gca,'XScale','log')
xlabel('Timestep')
ylabel('p(t_{D} \geq t)')
title('CDF');
legend('L', 'R', 'Location', 'southeast')

%% r traces (not for presentation)
figure; 
subplot(1,2,1); hold on; 
subplot(1,2,2); hold on; 
for q = 1:250
    rs = r{q}';
    if anyDecision(q)
        subplot(1,2,1)
        plot(rs(:, 1), 'r')
    	plot(rs(:, 2), 'b')
    else
        subplot(1,2,2)
        plot(rs(:, 1), 'r')
        plot(rs(:, 2), 'b')
    end

end

%% r traces + DT histogram
nBins1 = maxsteps + 1;  % decision time histogram
nBins2 = 20;  % undecided ending-r histogram

% plot parameters
set(0, 'defaultAxesFontSize', 12);
xCutoff = maxsteps;
Nplot = 100;  % number of r traces to plot
fracL = wins(1) / nReps;
fracR = wins(2) / nReps;
lineAlpha = 0.4;  % transparency of r traces
withText = false;  
textSize = 10;
colors = lines(2);
c1 = colors(1, :); 
c2 = colors(2, :);

% subplot positions
left = 0.1;
top = 0.9;
bottom1 = 0.70;
bottom2 = 0.15;
width = 0.85;
right = 0.95;
height=0.2;
histHeight = top - bottom1;
mainHeight = bottom1 - bottom2;
right1 = right - histHeight;

% create figure with size and colour
fig1 = figure('Color', 'w', 'Position', [800 100 750 600]);

% histograms of decision times
ax1 = axes('Position', [left bottom1 right1-left histHeight]);
ylabel('p(D)')
rtNonZero = rt;
rtNonZero(rtNonZero == 0) = nan; %maxsteps+1;
hist1 = histogram(rtNonZero(:, 1), 'Normalization', 'pdf'); % 'DisplayStyle', 'stairs');%, nBins1);
hold on;
hist2 = histogram(rtNonZero(:, 2), hist1.BinEdges, 'Normalization', 'pdf'); %, 'DisplayStyle', 'stairs');
%[rtCounts, rtCenters] = hist(rtNonZero, nBins1);
%totalCount = sum(cumsum(rtCounts), 2);
%xCutoff = find(totalCount > totalCount(end) - 1, 1);
%bar1 = bar(rtCenters, rtCounts / totalCount(end), 1);
%bar1(1).FaceColor = c1;
%bar1(2).FaceColor = c2;
hist1.FaceColor = c1;
hist2.FaceColor = c2;
hist1.LineStyle = 'none';
hist2.LineStyle = 'none';
%hist1.LineWidth = 2;
%hist2.LineWidth = 2;
%hist1.DisplayStyle = 'stairs';
%hist2.DisplayStyle = 'stairs';
xlim([0 xCutoff])
set(gca, 'XTickLabel', [],'XTick',[], 'YTick', [], 'box', 'off')
hAx = gca;
hAx.YRuler.Axle.LineStyle = 'none';
legend({'L', 'R'}, 'Location', 'northeast', 'FontSize',11)
legend('boxoff')
nDecided = sum(anyDecision);
nUndecided = nReps - nDecided;
%maxCount = max(max(hist1.BinCounts/nDecided), max(hist2.BinCounts/nDecided));
%maxCount = max(max(hist1.BinCounts/sum(hist1.BinCounts)), max(hist2.BinCounts/sum(hist2.BinCounts)));

fitType = 'Stable';
hist1Fit = fitdist(rtNonZero(:, 1), fitType); 
%hist1FitKernel = fitdist(rtNonZero(:, 1), 'Kernel'); 
%hist1FitIG = fitdist(rtNonZero(:, 1), 'InverseGaussian');
hist2Fit = fitdist(rtNonZero(:, 2), fitType);
%hist2FitKernel = fitdist(rtNonZero(:, 2), 'Kernel'); 
%hist2FitIG = fitdist(rtNonZero(:, 2), 'InverseGaussian');
plot(1:xCutoff, pdf(hist1Fit, 1:xCutoff) * fracL, 'LineWidth', 1.5, 'Color', c1, 'HandleVisibility','off')
plot(1:xCutoff, pdf(hist2Fit, 1:xCutoff) * fracR, 'LineWidth', 1.5, 'Color', c2, 'HandleVisibility','off')
plot(1:xCutoff, pdf(hist1FitIG, 1:xCutoff) * fracL, 'LineStyle', ':', 'LineWidth', 1.5, 'Color', c1, 'HandleVisibility','off')
plot(1:xCutoff, pdf(hist2FitIG, 1:xCutoff) * fracR, 'LineStyle', ':', 'LineWidth', 1.5, 'Color', c2, 'HandleVisibility','off')
%plot(1:xCutoff, pdf(hist1FitKernel, 1:xCutoff), 'LineStyle', ':', 'LineWidth', 2, 'Color', c1, 'HandleVisibility','off')
%plot(1:xCutoff, pdf(hist2FitKernel, 1:xCutoff), 'LineStyle', ':', 'LineWidth', 2, 'Color', c2, 'HandleVisibility','off')


% r over time in each trial
ax2 = axes('Position', [left bottom2 right1-left mainHeight]); 
hold on;
NplotMin = min(Nplot, size(r, 2));
sampleTrials = randsample(nReps, NplotMin);
rs_all = nan(NplotMin, 2, maxsteps);
for q = 1:NplotMin
    rs_all(q, :, 1:size(r{sampleTrials(q)}, 2)) = r{sampleTrials(q)};
end
for q = 1:NplotMin
    rs = squeeze(rs_all(q, :, :));
    plotSteps = size(rs, 2);
    plot1 = plot(1:plotSteps, rs(1, :));
    plot2 = plot(1:plotSteps, rs(2, :));
    plot1.Color(1:3) = c1;
    plot1.Color(4) = lineAlpha;
    plot1.LineWidth = 1;
    plot2.Color(1:3) = c2;
    plot2.Color(4) = lineAlpha;
    plot1.LineWidth = 1;
end
xlim([0 xCutoff])
xlabel('Timestep')
ylabel('r')
hAx = gca;
hAx.BoxStyle = 'full';
ylim([-0 thres])
set(gca, 'YTick', [0 0.5 thres])
hAx.YRuler.Axle.LineStyle = 'solid';
hAx.XRuler.Axle.LineStyle = 'solid';
box on;
if withText
    text(xCutoff*0.99, thres*0.95, sprintf('N_{plot} = %i', firstN), 'HorizontalAlignment', 'right', 'FontSize', textSize);
end

% histogram of r at end for undecided trials
ax3 = axes('Position', [right1 bottom2 histHeight mainHeight]); 
undecidedTrials = find(nansum(rtNonZero, 2) == 0);
undecidedEndRs = [];
for u = 1:length(undecidedTrials)
    undecidedEndRs = [undecidedEndRs; r{undecidedTrials(u)}(:, end)'];
end
%[undecidedRCounts, undecidedRCenters] = hist(undecidedEndRs, nBins2);
%totalCount = sum(cumsum(undecidedRCounts), 2);
hist3 = histogram(undecidedEndRs(:, 1), nBins2, 'Normalization', 'pdf', 'Orientation', 'horizontal'); % 'DisplayStyle', 'stairs', 
hold on;
hist4 = histogram(undecidedEndRs(:, 2), hist3.BinEdges, 'Normalization', 'pdf', 'Orientation', 'horizontal'); % 'DisplayStyle', 'stairs', 
%bar2 = barh(undecidedRCenters, undecidedRCounts / totalCount(end), 1);
%bar2(1).FaceColor = c1;
%bar2(2).FaceColor = c2;
%bar2(1).EdgeColor = c1;
%bar2(2).EdgeColor = c2;
set(gca, 'XTickLabel', [],'XTick',[], 'YTick', [], 'box', 'off')
ylim([-0 thres])
hAx = gca;
hist3.FaceColor = c1;
hist4.FaceColor = c2;
hist3.LineStyle = 'none';
hist4.LineStyle = 'none';
%hist3.LineWidth = 2;
%hist4.LineWidth = 2;
%hAx.YRuler.Axle.LineStyle = 'none';
hist3Fit = fitdist(undecidedEndRs(:, 1), 'Exponential');
hist4Fit = fitdist(undecidedEndRs(:, 2), 'Exponential');
plot(pdf(hist3Fit, 0:0.01:thres), 0:0.01:thres, 'LineWidth', 1.5, 'Color', c1, 'HandleVisibility','off')
plot(pdf(hist4Fit, 0:0.01:thres), 0:0.01:thres, 'LineWidth', 1.5, 'Color', c2, 'HandleVisibility','off')
hAx.XAxis.Visible = false;
maxUndecidedCount = max(max(hist3.BinCounts), max(hist4.BinCounts));
if withText
    text(maxUndecidedCount * 0.05, thres*0.95, sprintf('N_{?} = %i', nUndecided), 'HorizontalAlignment', 'left', 'FontSize', textSize);
    text(maxUndecidedCount * 0.05, thres*(1 + 0.7*histHeight/mainHeight), sprintf('N = %i', nReps), 'HorizontalAlignment', 'left', 'FontSize', textSize);
end
saveas(fig1, 'rs-with-histograms.png')


%% sequential decision analysis
% decision = zeros(seqTrials, nNets, nSeq, 2);
figure;
q = 1;  % which repeat to plot
decisions = squeeze(decision(q, :, :, 1));  
subplot(2,1,1); hold on;
plot(1:seqTrials, trueSeq{q}, '*k');
plot(1:seqTrials, decisions(1, :), '.b')
plot(1:seqTrials, decisions(2, :), '.r')
ylim([1 2])
subplot(2,1,2);
plot(1:seqTrials, lams{q}(:, :, 1))