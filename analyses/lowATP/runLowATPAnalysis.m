% Effect of reduced intracellular ATP (8 -> 4 -> 2 mM) on the Ca transient of the
% T-World endocardial cell, with free ADP estimated from the creatine kinase (CK)
% equilibrium. See README.md in this folder for the write-up.
%
% Produces figures/*.png and results.csv in this folder.
clear

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
figDir = fullfile(thisDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

beats = 300;  % 5 min at 1 Hz; per-beat summaries are kept for all beats
bcl = 1000;
X0 = getStartingState('TW_endo');

%% Energetic scenarios
% PCr buffers ATP via CK, so [ATP] falls only once PCr is largely depleted. The
% residual PCr fraction and pH at each ATP level are assumptions (see README).
ATP_ref = 8;
[ADP_ref, ref] = getCKEquilibriumADP(ATP_ref, 1, 7.1);

sc = struct('name', {}, 'label', {}, 'ATP', {}, 'fPCr', {}, 'pH', {}, 'ADP', {}, 'IKATP', {});
sc(end+1) = struct('name', 'ATP8', 'label', 'ATP 8 mM (reference)', 'ATP', 8, 'fPCr', 1, 'pH', 7.1, 'ADP', ADP_ref, 'IKATP', 1);
lowATP = [4, 2];
lowPCr = [0.10, 0.05];
lowpH_acid = [6.4, 6.2];
for k = 1:2
    a = lowATP(k);
    sc(end+1) = struct('name', sprintf('ATP%d_acid', a), 'label', sprintf('ATP %d mM, CK ischaemic (pH %.1f)', a, lowpH_acid(k)), ...
        'ATP', a, 'fPCr', lowPCr(k), 'pH', lowpH_acid(k), 'ADP', getCKEquilibriumADP(a, lowPCr(k), lowpH_acid(k)), 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_neutral', a), 'label', sprintf('ATP %d mM, CK non-acidotic (pH 7.0)', a), ...
        'ATP', a, 'fPCr', lowPCr(k), 'pH', 7.0, 'ADP', getCKEquilibriumADP(a, lowPCr(k), 7.0), 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_ADPref', a), 'label', sprintf('ATP %d mM, ADP at reference', a), ...
        'ATP', a, 'fPCr', NaN, 'pH', NaN, 'ADP', ADP_ref, 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_acid_noKATP', a), 'label', sprintf('ATP %d mM, CK ischaemic, I_{KATP} blocked', a), ...
        'ATP', a, 'fPCr', lowPCr(k), 'pH', lowpH_acid(k), 'ADP', getCKEquilibriumADP(a, lowPCr(k), lowpH_acid(k)), 'IKATP', 0); %#ok<SAGROW>
end

%% Simulations
nSc = numel(sc);
res = struct();
for s = 1:nSc
    param = [];
    param.model = @model_TWorld;
    param.bcl = bcl;
    param.ATP_i = sc(s).ATP;
    param.ADP_i = sc(s).ADP;
    param.IKATP_Multiplier = sc(s).IKATP;
    fprintf('Simulating %s (ADP = %.1f uM)\n', sc(s).label, 1000 * sc(s).ADP);

    [time, X] = modelRunner(X0, [], param, beats, 0);

    % per-beat summaries (Ca transient amplitude and [Na]i drift)
    perBeat = zeros(beats, 3);
    for b = 1:beats
        cai = X{b}(:, 8);
        perBeat(b, :) = [max(cai) - min(cai), min(cai), X{b}(end, 4)];
    end

    currents = getCurrentsStructure(time(end), X(end), param, 0);
    m = getBeatMetrics(currents);
    m.perBeat = perBeat;
    m.currents = currents;
    m.factors = [currents.fATP_SERCA(1), currents.fATP_NaK(1), currents.fATP_PCa(1), currents.fATP_KATP(1)];
    res.(sc(s).name) = m;
end

%% Results table
names = {sc.name}';
T = table(names, {sc.label}', [sc.ATP]', 1000 * [sc.ADP]', [sc.fPCr]', [sc.pH]', [sc.IKATP]', ...
    'VariableNames', {'scenario', 'label', 'ATP_mM', 'ADP_uM', 'fPCr', 'pH', 'IKATP_on'});
metricNames = {'CaT_amp_nM', 'CaT_dia_nM', 'CaT_peak_nM', 'CaT_ttp_ms', 'CaT_tau_ms', 'CaTD50_ms', 'CaTD90_ms', ...
    'APD90_ms', 'V_rest_mV', 'CaSR_max_mM', 'Nai_mM', 'Ta_peak_kPa', 'IKATP_peak'};
for k = 1:numel(metricNames)
    T.(metricNames{k}) = cellfun(@(n) res.(n).(metricNames{k}), names);
end
T.fSERCA = cellfun(@(n) res.(n).factors(1), names);
T.fNaK = cellfun(@(n) res.(n).factors(2), names);
T.fPCa = cellfun(@(n) res.(n).factors(3), names);
T.fKATP_open = cellfun(@(n) res.(n).factors(4), names);
writetable(T, fullfile(thisDir, 'results.csv'));
disp(T(:, [1 3 4 8 9 10 12 13 14 16 17 19 20 21 22 23 24]));

%% Figures
cols = lines(7);
setStyle = @(ax) set(ax, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

% 1) Main comparison: Ca transient, AP, tension for the central (ischaemic CK) scenarios
main = {'ATP8', 'ATP4_acid', 'ATP2_acid'};
mainLabels = {'8 mM (ref.)', '4 mM', '2 mM'};
f = figure('Position', [100 100 1100 330], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
panels = {'Ca_i', 1e6, '[Ca^{2+}]_i (nM)'; 'V', 1, 'V_m (mV)'; 'Land_Ta', 1, 'Active tension (kPa)'};
for p = 1:3
    ax = nexttile(tl); hold(ax, 'on'); setStyle(ax);
    for k = 1:numel(main)
        c = res.(main{k}).currents;
        plot(ax, c.time - c.time(1), c.(panels{p, 1}) * panels{p, 2}, 'LineWidth', 1.5, 'Color', cols(k, :));
    end
    xlim(ax, [0 600]); xlabel(ax, 'Time (ms)'); ylabel(ax, panels{p, 3});
end
legend(nexttile(tl, 1), mainLabels, 'Location', 'northeast', 'Box', 'off');
title(tl, 'ATP_i reduction with CK-equilibrium ADP (ischaemic pH), beat 300 at 1 Hz');
exportgraphics(f, fullfile(figDir, 'fig1_CaT_AP_tension.png'), 'Resolution', 200);

% 2) Normalised Ca transients (shape / kinetics)
f = figure('Position', [100 100 520 380], 'Color', 'w'); ax = axes(f); hold(ax, 'on'); setStyle(ax);
for k = 1:numel(main)
    c = res.(main{k}).currents; cai = c.Ca_i;
    plot(ax, c.time - c.time(1), (cai - min(cai)) / (max(cai) - min(cai)), 'LineWidth', 1.5, 'Color', cols(k, :));
end
xlim(ax, [0 800]); xlabel(ax, 'Time (ms)'); ylabel(ax, 'Normalised \Delta[Ca^{2+}]_i');
legend(ax, mainLabels, 'Box', 'off'); title(ax, 'Ca transient shape');
exportgraphics(f, fullfile(figDir, 'fig2_CaT_normalised.png'), 'Resolution', 200);

% 3) Mechanism dissection at 4 and 2 mM
f = figure('Position', [100 100 1000 380], 'Color', 'w');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
for k = 1:2
    a = lowATP(k);
    ax = nexttile(tl); hold(ax, 'on'); setStyle(ax);
    set_ = {'ATP8', sprintf('ATP%d_acid', a), sprintf('ATP%d_neutral', a), sprintf('ATP%d_ADPref', a), sprintf('ATP%d_acid_noKATP', a)};
    lab = {'8 mM (ref.)', sprintf('%d mM, ADP %.0f \\muM (pH %.1f)', a, 1000 * getCKEquilibriumADP(a, lowPCr(k), lowpH_acid(k)), lowpH_acid(k)), ...
        sprintf('%d mM, ADP %.0f \\muM (pH 7.0)', a, 1000 * getCKEquilibriumADP(a, lowPCr(k), 7.0)), ...
        sprintf('%d mM, ADP %.0f \\muM (ref.)', a, 1000 * ADP_ref), sprintf('%d mM, pH %.1f, I_{KATP} blocked', a, lowpH_acid(k))};
    styles = {'-', '-', '--', ':', '-.'};
    colIdx = [1, k + 1, k + 1, k + 1, 5];
    for j = 1:numel(set_)
        c = res.(set_{j}).currents;
        plot(ax, c.time - c.time(1), c.Ca_i * 1e6, styles{j}, 'LineWidth', 1.5, 'Color', cols(colIdx(j), :));
    end
    xlim(ax, [0 700]); xlabel(ax, 'Time (ms)'); ylabel(ax, '[Ca^{2+}]_i (nM)');
    title(ax, sprintf('[ATP]_i = %d mM', a)); legend(ax, lab, 'Box', 'off', 'FontSize', 8);
end
exportgraphics(f, fullfile(figDir, 'fig3_mechanism_dissection.png'), 'Resolution', 200);

% 4) Beat-to-beat evolution of CaT amplitude and [Na]i
f = figure('Position', [100 100 1000 330], 'Color', 'w');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax1 = nexttile(tl); hold(ax1, 'on'); setStyle(ax1);
ax2 = nexttile(tl); hold(ax2, 'on'); setStyle(ax2);
for k = 1:numel(main)
    pb = res.(main{k}).perBeat;
    plot(ax1, 1:beats, pb(:, 1) * 1e6, 'LineWidth', 1.5, 'Color', cols(k, :));
    plot(ax2, 1:beats, pb(:, 3), 'LineWidth', 1.5, 'Color', cols(k, :));
end
xlabel(ax1, 'Beat'); ylabel(ax1, 'CaT amplitude (nM)'); xlabel(ax2, 'Beat'); ylabel(ax2, '[Na^+]_i (mM)');
legend(ax1, mainLabels, 'Box', 'off', 'Location', 'southeast');
ylim(ax1, [270 340]);
exportgraphics(f, fullfile(figDir, 'fig4_beat_evolution.png'), 'Resolution', 200);

% 5) Nucleotide dependence of the model components
ATPgrid = logspace(log10(0.5), log10(10), 200);
ADPlevels = [ADP_ref, 0.1, 0.4];
f = figure('Position', [100 100 1100 330], 'Color', 'w');
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact');
titles = {'SERCA', 'Na^+/K^+ pump', 'PMCA (I_{pCa})', 'K_{ATP} open fraction (endo)'};
for p = 1:4
    ax = nexttile(tl); hold(ax, 'on'); setStyle(ax); set(ax, 'XScale', 'log');
    for j = 1:numel(ADPlevels)
        y = zeros(size(ATPgrid));
        for i = 1:numel(ATPgrid)
            y(i) = getFactor(p, ATPgrid(i), ADPlevels(j), ATP_ref, ADP_ref);
        end
        plot(ax, ATPgrid, y, 'LineWidth', 1.5, 'Color', cols(j, :));
    end
    xline(ax, [2 4 8], ':', 'Color', [0.5 0.5 0.5]);
    xlabel(ax, '[ATP]_i (mM)'); title(ax, titles{p});
    if p < 4, ylabel(ax, 'Rate relative to reference'); else, ylabel(ax, 'Excess open fraction'); end
end
legend(nexttile(tl, 1), arrayfun(@(d) sprintf('ADP %.0f \\muM', 1000 * d), ADPlevels, 'UniformOutput', false), 'Box', 'off', 'Location', 'southeast');
exportgraphics(f, fullfile(figDir, 'fig5_nucleotide_factors.png'), 'Resolution', 200);

% 6) CK-equilibrium free ADP vs residual PCr and pH
fPCrGrid = linspace(0.02, 0.5, 100);
pHGrid = linspace(6.0, 7.2, 100);
f = figure('Position', [100 100 900 360], 'Color', 'w');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
for k = 1:2
    a = lowATP(k);
    Z = zeros(numel(pHGrid), numel(fPCrGrid));
    for i = 1:numel(pHGrid)
        for j = 1:numel(fPCrGrid)
            Z(i, j) = 1000 * getCKEquilibriumADP(a, fPCrGrid(j), pHGrid(i));
        end
    end
    ax = nexttile(tl); setStyle(ax);
    contourf(ax, 100 * fPCrGrid, pHGrid, log10(Z), 20, 'LineStyle', 'none'); hold(ax, 'on');
    [C, h] = contour(ax, 100 * fPCrGrid, pHGrid, Z, [25 50 100 200 400 800 1600], 'k');
    clabel(C, h, 'FontSize', 8);
    plot(ax, 100 * lowPCr(k), lowpH_acid(k), 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    plot(ax, 100 * lowPCr(k), 7.0, 'ws', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
    xlabel(ax, 'Residual PCr (% of normoxic)'); ylabel(ax, 'pH_i');
    title(ax, sprintf('Free [ADP] (\\muM) at [ATP]_i = %d mM', a));
    colormap(ax, parula);
    caxis(ax, log10([10 3000]));
    cb = colorbar(ax);
    cb.Ticks = log10([10 30 100 300 1000 3000]);
    cb.TickLabels = {'10', '30', '100', '300', '1000', '3000'};
    cb.Label.String = 'Free [ADP] (\muM)';
end
exportgraphics(f, fullfile(figDir, 'fig6_CK_ADP_map.png'), 'Resolution', 200);

save(fullfile(thisDir, 'results_summary.mat'), 'T', 'sc', 'ref');

%% Local functions
function m = getBeatMetrics(c)
t = c.time - c.time(1);
cai = c.Ca_i * 1e6;
[caPeak, iPeak] = max(cai);
caDia = min(cai);
amp = caPeak - caDia;
m.CaT_amp_nM = amp;
m.CaT_dia_nM = caDia;
m.CaT_peak_nM = caPeak;
m.CaT_ttp_ms = t(iPeak);
% single-exponential decay fitted from 10 % decay (excluding the rounded peak) on
iStart = iPeak - 1 + find(cai(iPeak:end) <= caDia + 0.9 * amp, 1);
tFit = t(iStart:end) - t(iStart);
cFit = cai(iStart:end);
sse = @(p) sum((p(1) * exp(-tFit / p(2)) + p(3) - cFit).^2);
p = fminsearch(sse, [cFit(1) - cFit(end), 100, cFit(end)], optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4));
m.CaT_tau_ms = p(2);
m.CaTD50_ms = crossingTime(t, cai, iPeak, caPeak - 0.5 * amp);
m.CaTD90_ms = crossingTime(t, cai, iPeak, caPeak - 0.9 * amp);
v = c.V;
[vPeak, iv] = max(v);
vThresh = vPeak - 0.9 * (vPeak - v(1));
m.APD90_ms = crossingTime(t, v, iv, vThresh) - interp1(v(1:iv), t(1:iv), vThresh);
m.V_rest_mV = v(1);
m.CaSR_max_mM = max(c.Ca_SR);
m.Nai_mM = c.Na_i(end);
m.Ta_peak_kPa = max(c.Land_Ta);
m.IKATP_peak = max(c.IKATP);
end

function tc = crossingTime(t, y, iPeak, level)
% time at which y falls below level after its peak (linear interpolation)
i = iPeak - 1 + find(y(iPeak:end) <= level, 1);
tc = interp1(y([i - 1, i]), t([i - 1, i]), level);
end

function y = getFactor(component, ATP, ADP, ATP_ref, ADP_ref)
% Normalised nucleotide factors as implemented in model_TWorld (endocardial cell).
g = @(a, d) [1 / (0.01/a * (1 + d/0.14) + 1 + d/5.1), ...
    1 / (1 + 0.008/a * (1 + d/0.1)), ...
    1 / (1 + 0.012/a * (1 + d/1.0)) + 1 / (1 + 0.23/a), ...
    1 / (1 + (a * 1000 / (0.32 * (35.8 + 17.9 * (d * 1000)^0.256)))^(1.3 + 4.44 * exp(-0.09 * d * 1000)))];
gi = g(ATP, ADP);
gr = g(ATP_ref, ADP_ref);
if component < 4
    y = gi(component) / gr(component);
else
    y = max(0, (gi(4) - gr(4)) / (1 - gr(4)));
end
end
