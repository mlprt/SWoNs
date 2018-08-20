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
lineAlpha = 0.35;
colors = lines(2);
c1 = colors(1, :);
c2 = colors(2, :);

% subplot positions
left = 0.1;
top = 0.9;
bottom1 = 0.75;
bottom2 = 0.15;
width=0.85;
right = 0.95;
height=0.2;
histHeight = top - bottom1;
mainHeight = bottom1 - bottom2;
right1 = 0.95 - histHeight;

fig1 = figure('Color', 'w', 'Position', [800 100 750 600]);

% histograms of decision times
ax1 = axes('Position', [left bottom1 right1-left histHeight]);
ylabel('p(D)')
rtNonZero = rt;
rtNonZero(rtNonZero == 0) = nan; %maxsteps+1;
hist1 = histogram(rtNonZero(:, 1), 'Normalization', 'probability'); % 'DisplayStyle', 'stairs');%, nBins1);
hold on;
hist2 = histogram(rtNonZero(:, 2), hist1.NumBins, 'Normalization', 'probability'); %, 'DisplayStyle', 'stairs');
%[rtCounts, rtCenters] = hist(rtNonZero, nBins1);
xCutoff = maxsteps;
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
legend('L', 'R', 'Location', 'northeast')
legend('boxoff')

% r over time in each trial
ax2 = axes('Position', [left bottom2 right1-left mainHeight]); 
hold on;
rs = nan(size(r, 2), 2, maxsteps);
for q = 1:size(r, 2)
    rs(q, :, 1:size(r{q}, 2)) = r{q};
end
plot1 = plot(1:maxsteps, squeeze(rs(:, 1, :)), 'Color', [c1 lineAlpha]);
plot2 = plot(1:maxsteps, squeeze(rs(:, 2, :)), 'Color', [c2 lineAlpha]);
%plot1.Color(:, :) = ;
%plot2.Color(:, :) = [c2 lineAlpha];
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

% histogram of r at end for undecided trials
ax3 = axes('Position', [right1 bottom2 histHeight mainHeight]); 
undecidedTrials = find(nansum(rtNonZero, 2) == 0);
undecidedEndRs = [];
for u = 1:length(undecidedTrials)
    undecidedEndRs = [undecidedEndRs; r{undecidedTrials(u)}(:, end)'];
end
%[undecidedRCounts, undecidedRCenters] = hist(undecidedEndRs, nBins2);
%totalCount = sum(cumsum(undecidedRCounts), 2);
hist3 = histogram(undecidedEndRs(:, 1), nBins2, 'Orientation', 'horizontal'); % 'DisplayStyle', 'stairs', 
hold on;
hist4 = histogram(undecidedEndRs(:, 2), hist3.BinEdges, 'Orientation', 'horizontal'); % 'DisplayStyle', 'stairs', 
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
hAx.XAxis.Visible = false;
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