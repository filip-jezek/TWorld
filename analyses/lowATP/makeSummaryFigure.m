% Summary figure for SUMMARY.md: Ca transient and Land active tension at 8, 4 and 2 mM
% ATP (lower-bound scenarios, lowest published ischaemic ADP), and their % changes.
% Traces are re-simulated (300 beats, last beat shown); bar values and the range
% (I_KATP g0 x0.5 to highest published ADP) are read from results.csv, so run
% runLowATPAnalysis first.
clear

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
figDir = fullfile(thisDir, 'figures');
T = readtable(fullfile(thisDir, 'results.csv'), 'TextType', 'string');

ATP = [8, 4, 2];
ADP = [getCKEquilibriumADP(8), 0.182, 0.109]; % mM; reference CK value, lowest published ischaemic
beats = 300;
X0 = getStartingState('TW_endo');

cacheFile = fullfile(tempdir, 'lowATP_summary_traces.mat'); % avoids re-simulating when only restyling
if isfile(cacheFile)
    load(cacheFile, 'traces');
else
    traces = cell(1, numel(ATP));
    for k = 1:numel(ATP)
        param = [];
        param.model = @model_TWorld;
        param.bcl = 1000;
        param.ATP_i = ATP(k);
        param.ADP_i = ADP(k);
        fprintf('Simulating ATP %g mM\n', ATP(k));
        [time, X] = modelRunner(X0, [], param, beats, beats - 1);
        traces{k} = getCurrentsStructure(time, X, param, 0);
    end
    save(cacheFile, 'traces');
end

% % changes vs 8 mM: lower bound, and the range from g0 x0.5 (smallest) to highest ADP (largest)
row = @(name) T(T.scenario == name, :);
ref = row("ATP8");
pct = @(r, v) 100 * (r.(v) / ref.(v) - 1);
groups = {'CaT_amp_nM', 'Ta_peak_kPa'};
lb = zeros(2, 2); lo = lb; hi = lb; % rows: 4, 2 mM; cols: CaT, Ta
for i = 1:2
    a = ATP(i + 1);
    for g = 1:2
        lb(i, g) = pct(row(sprintf("ATP%d_ADPmin", a)), groups{g});
        lo(i, g) = pct(row(sprintf("KATP_ATP%d_ADPmin_g050_Km100", a)), groups{g});
        hi(i, g) = pct(row(sprintf("ATP%d_ADPmax", a)), groups{g});
    end
end

cols = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125];
labels = arrayfun(@(a, d) sprintf('%g mM ATP (ADP %.0f \\muM)', a, 1000 * d), ATP, ADP, 'UniformOutput', false);
f = figure('Position', [100 100 1150 340], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
panels = {'Ca_i', 1e6, '[Ca^{2+}]_i (nM)', 'Ca^{2+} transient'; 'Land_Ta', 1, 'Active tension (kPa)', 'Land active tension'};
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
    for k = 1:numel(ATP)
        c = traces{k};
        plot(ax, c.time - c.time(1), c.(panels{p, 1}) * panels{p, 2}, 'LineWidth', 1.6, 'Color', cols(k, :));
    end
    xlim(ax, [0 600]); xlabel(ax, 'Time (ms)'); ylabel(ax, panels{p, 3}); title(ax, panels{p, 4});
end
legend(nexttile(tl, 1), labels, 'Box', 'off', 'FontSize', 8);

ax = nexttile(tl); hold(ax, 'on'); set(ax, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
b = bar(ax, 1:2, lb', 'grouped');
for i = 1:2
    b(i).FaceColor = cols(i + 1, :);
    x = b(i).XEndPoints;
    errorbar(ax, x, lb(i, :), lb(i, :) - hi(i, :), lo(i, :) - lb(i, :), 'k', 'LineStyle', 'none', 'LineWidth', 1);
    for g = 1:2
        text(ax, x(g), hi(i, g) - 1, sprintf('%.1f%%', lb(i, g)), 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
end
set(ax, 'XTick', 1:2, 'XTickLabel', {'CaT amplitude', 'Peak tension'});
ylabel(ax, 'Change vs 8 mM ATP (%)'); yline(ax, 0, 'k');
xlim(ax, [0.5 2.9]); ylim(ax, [5 * floor((min(hi(:)) - 5) / 5), 0]);
legend(ax, b, {'4 mM', '2 mM'}, 'Box', 'off', 'Location', 'southwest');
title(ax, {'Lower bound (bars, labelled)', 'whiskers: I_{KATP} g_0 \times0.5 to highest published ADP'}, 'FontSize', 9);
title(tl, 'Endocardial T-World, 1 Hz, beat 300; ATP and free ADP clamped');
exportgraphics(f, fullfile(figDir, 'fig0_summary.png'), 'Resolution', 200);
