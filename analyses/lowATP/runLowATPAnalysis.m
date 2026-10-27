% Effect of reduced intracellular ATP (8 -> 4 -> 2 mM) on the Ca transient of the
% T-World endocardial cell, aimed at the minimal (lower-bound) effect: free ADP at
% reduced ATP is taken from published ischaemic data (PCr is exhausted, so the creatine
% kinase equilibrium is ill-conditioned), and the normoxic reference ADP from the human
% CK equilibrium. See README.md in this folder for the write-up.
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
% Reference: human normoxic pools in CK equilibrium (Weiss et al. 2005), free ADP ~94 uM.
% A high reference ADP minimises the computed change, so it is the conservative choice.
ATP_ref = 8;
[ADP_ref, ref] = getCKEquilibriumADP(ATP_ref);

% Published free ADP during global no-flow ischaemia at matching ATP fractions (50 %, 25 %):
%  min: guinea pig, Ferrero et al. 2023 Fig. 1 (Terkildsen 2007 profile), interpolated to
%       ATP 50 % (182 uM) and 25 % (109 uM) of control;
%  max: rat, Bak & Ingwall 1994 Table II, CK-calculated (284 uM at 54 % ATP, 767 uM at 31 %).
lowATP = [4, 2];
ADPmin = [0.182, 0.109];
ADPmax = [0.284, 0.767];

sc = struct('name', {}, 'label', {}, 'ATP', {}, 'fPCr', {}, 'pH', {}, 'ADP', {}, 'IKATP', {});
sc(end+1) = struct('name', 'ATP8', 'label', 'ATP 8 mM (reference)', 'ATP', 8, 'fPCr', 1, 'pH', ref.pH, 'ADP', ADP_ref, 'IKATP', 1);
for k = 1:2
    a = lowATP(k);
    sc(end+1) = struct('name', sprintf('ATP%d_ADPmin', a), 'label', sprintf('ATP %d mM, ADP %.0f uM (lowest published)', a, 1000 * ADPmin(k)), ...
        'ATP', a, 'fPCr', NaN, 'pH', NaN, 'ADP', ADPmin(k), 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_ADPmax', a), 'label', sprintf('ATP %d mM, ADP %.0f uM (highest published)', a, 1000 * ADPmax(k)), ...
        'ATP', a, 'fPCr', NaN, 'pH', NaN, 'ADP', ADPmax(k), 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_ADPref', a), 'label', sprintf('ATP %d mM, ADP unchanged (%.0f uM)', a, 1000 * ADP_ref), ...
        'ATP', a, 'fPCr', NaN, 'pH', NaN, 'ADP', ADP_ref, 'IKATP', 1); %#ok<SAGROW>
    sc(end+1) = struct('name', sprintf('ATP%d_ADPmin_noKATP', a), 'label', sprintf('ATP %d mM, ADP %.0f uM, I_{KATP} blocked', a, 1000 * ADPmin(k)), ...
        'ATP', a, 'fPCr', NaN, 'pH', NaN, 'ADP', ADPmin(k), 'IKATP', 0); %#ok<SAGROW>
end
[sc.KmMult] = deal(1);   % K_ATP Km (alpha) multiplier
[sc.ca50Mult] = deal(1); % myofilament (TnC) ca50 multiplier
[sc.beats] = deal(beats);

%% Sensitivity scenarios for the lower bound (CaT amplitude stable after ~50 beats)
sweepBeats = 150;
% (a) K_ATP strength: g0 x {0.5, 1, 2} and Km x {1 (endo alpha 0.32), 3.125 (alpha 1, epi/Ferrero 1996)}
g0Mult = [0.5, 1, 2];
KmMult = [1, 3.125];
for k = 1:2
    a = lowATP(k);
    ADPs = [ADP_ref, ADPmin(k)];
    ADPtags = {'ADPref', 'ADPmin'};
    for d = 1:2
        for g = g0Mult
            for km = KmMult
                if g == 1 && km == 1, continue; end % already simulated above
                sc(end+1) = struct('name', sprintf('KATP_ATP%d_%s_g%03d_Km%03d', a, ADPtags{d}, round(100 * g), round(100 * km)), ...
                    'label', sprintf('ATP %d mM, %s, g0 x%.2g, Km x%.3g', a, ADPtags{d}, g, km), 'ATP', a, 'fPCr', NaN, 'pH', NaN, ...
                    'ADP', ADPs(d), 'IKATP', g, 'KmMult', km, 'ca50Mult', 1, 'beats', sweepBeats); %#ok<SAGROW>
            end
        end
    end
end
% (b) myofilament Ca2+ sensitivity: TnC ca50 x {0.9, 1.25, 1.5, 2} (dpCa50 +0.05 ... -0.30)
ca50Mult = [0.9, 1.25, 1.5, 2];
baseSc = {'ATP8', 'ATP4_ADPmin', 'ATP2_ADPmin'};
for b = 1:numel(baseSc)
    s0 = sc(strcmp({sc.name}, baseSc{b}));
    for cm = ca50Mult
        s1 = s0;
        s1.name = sprintf('ca50_%s_x%03d', baseSc{b}, round(100 * cm));
        s1.label = sprintf('%s, ca50 x%.2g', s0.label, cm);
        s1.ca50Mult = cm;
        s1.beats = sweepBeats;
        sc(end+1) = s1; %#ok<SAGROW>
    end
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
    param.KATP_Km_Multiplier = sc(s).KmMult;
    param.ca50_Multiplier = sc(s).ca50Mult;
    fprintf('Simulating %s (ADP = %.1f uM)\n', sc(s).label, 1000 * sc(s).ADP);

    [time, X] = modelRunner(X0, [], param, sc(s).beats, 0);

    % per-beat summaries (Ca transient amplitude and [Na]i drift)
    perBeat = zeros(sc(s).beats, 3);
    for b = 1:sc(s).beats
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
T = table(names, {sc.label}', [sc.ATP]', 1000 * [sc.ADP]', [sc.fPCr]', [sc.pH]', [sc.IKATP]', [sc.KmMult]', [sc.ca50Mult]', [sc.beats]', ...
    'VariableNames', {'scenario', 'label', 'ATP_mM', 'ADP_uM', 'fPCr', 'pH', 'IKATP_g0_mult', 'KATP_Km_mult', 'ca50_mult', 'beats'});
metricNames = {'CaT_amp_nM', 'CaT_dia_nM', 'CaT_peak_nM', 'CaT_ttp_ms', 'CaT_tau_ms', 'CaTD50_ms', 'CaTD90_ms', ...
    'APD90_ms', 'V_rest_mV', 'CaSR_max_mM', 'Nai_mM', 'Ta_peak_kPa', 'IKATP_peak'};
for k = 1:numel(metricNames)
    T.(metricNames{k}) = cellfun(@(n) res.(n).(metricNames{k}), names);
end
T.fSERCA = cellfun(@(n) res.(n).factors(1), names);
T.fNaK = cellfun(@(n) res.(n).factors(2), names);
T.fPCa = cellfun(@(n) res.(n).factors(3), names);
T.fKATP_open = cellfun(@(n) res.(n).factors(4), names);
T.dCaT_amp_pct = 100 * (T.CaT_amp_nM / T.CaT_amp_nM(1) - 1); % relative to ATP 8 mM reference
T.dAPD90_pct = 100 * (T.APD90_ms / T.APD90_ms(1) - 1);
writetable(T, fullfile(thisDir, 'results.csv'));
disp(T(:, {'scenario', 'ADP_uM', 'IKATP_g0_mult', 'KATP_Km_mult', 'ca50_mult', 'CaT_amp_nM', 'dCaT_amp_pct', 'CaT_tau_ms', 'CaTD90_ms', 'APD90_ms', 'dAPD90_pct'}));

%% Figures
cols = lines(7);
setStyle = @(ax) set(ax, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

% 1) Main comparison: Ca transient, AP, tension for the lower-bound (lowest published ADP) scenarios
main = {'ATP8', 'ATP4_ADPmin', 'ATP2_ADPmin'};
mainLabels = {sprintf('8 mM, ADP %.0f \\muM (ref.)', 1000 * ADP_ref), sprintf('4 mM, ADP %.0f \\muM', 1000 * ADPmin(1)), sprintf('2 mM, ADP %.0f \\muM', 1000 * ADPmin(2))};
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
title(tl, 'ATP_i reduction with lowest published ischaemic ADP (lower bound), beat 300 at 1 Hz');
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
    set_ = {'ATP8', sprintf('ATP%d_ADPmin', a), sprintf('ATP%d_ADPmax', a), sprintf('ATP%d_ADPref', a), sprintf('ATP%d_ADPmin_noKATP', a)};
    lab = {'8 mM (ref.)', sprintf('%d mM, ADP %.0f \\muM (lowest published)', a, 1000 * ADPmin(k)), ...
        sprintf('%d mM, ADP %.0f \\muM (highest published)', a, 1000 * ADPmax(k)), ...
        sprintf('%d mM, ADP unchanged (%.0f \\muM)', a, 1000 * ADP_ref), sprintf('%d mM, ADP %.0f \\muM, I_{KATP} blocked', a, 1000 * ADPmin(k))};
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
exportgraphics(f, fullfile(figDir, 'fig4_beat_evolution.png'), 'Resolution', 200);

% 5) Nucleotide dependence of the model components
ATPgrid = logspace(log10(0.5), log10(10), 200);
ADPlevels = [ADP_ref, ADPmin(1), ADPmax(2)];
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

% 6) Published free ADP vs ATP during global ischaemia, and the values adopted here
% guinea pig (Ferrero et al. 2023 Fig. 1, digitised, 0-30 min) and rat (Bak & Ingwall 1994 Table II)
gpATP = [10.0, 8.08, 6.28, 4.60, 3.13, 1.90, 1.04] / 10.0;
gpADP = [15, 328, 241, 164, 122, 96, 80];
ratATP = [10.8, 5.8, 3.4] / 10.8;
ratADP = [39, 284, 767];
f = figure('Position', [100 100 560 400], 'Color', 'w'); ax = axes(f); hold(ax, 'on'); setStyle(ax); set(ax, 'YScale', 'log');
plot(ax, 100 * gpATP, gpADP, '-o', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(ax, 100 * ratATP, ratADP, '-s', 'LineWidth', 1.2, 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.2 0.2 0.2]);
plot(ax, 100, 1000 * ADP_ref, 'p', 'MarkerSize', 12, 'MarkerFaceColor', cols(1, :), 'Color', cols(1, :));
plot(ax, 100 * lowATP / ATP_ref, 1000 * ADPmin, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cols(2, :), 'Color', cols(2, :));
plot(ax, 100 * lowATP / ATP_ref, 1000 * ADPmax, '^', 'MarkerSize', 9, 'MarkerFaceColor', cols(3, :), 'Color', cols(3, :));
set(ax, 'XDir', 'reverse'); xlim(ax, [0 105]);
xlabel(ax, '[ATP]_i (% of control)'); ylabel(ax, 'Free [ADP]_i (\muM)');
legend(ax, {'Guinea pig (Ferrero 2023 / Terkildsen 2007)', 'Rat (Bak & Ingwall 1994)', 'Reference, human CK (Weiss 2005)', ...
    'Adopted: lowest published', 'Adopted: highest published'}, 'Box', 'off', 'Location', 'south', 'FontSize', 8);
title(ax, 'Free ADP during global ischaemia');
exportgraphics(f, fullfile(figDir, 'fig6_ADP_vs_ATP_literature.png'), 'Resolution', 200);

% 7) Lower bound vs K_ATP strength (g0, Km)
getRow = @(n) T(strcmp(T.scenario, n), :);
f = figure('Position', [100 100 1000 380], 'Color', 'w');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
for k = 1:2
    a = lowATP(k);
    ax = nexttile(tl); hold(ax, 'on'); setStyle(ax); set(ax, 'XScale', 'log');
    lab = {};
    j = 0;
    for d = {'ADPref', 'ADPmin'}
        for km = KmMult
            j = j + 1;
            y = zeros(size(g0Mult));
            for gi = 1:numel(g0Mult)
                if g0Mult(gi) == 1 && km == 1
                    n = sprintf('ATP%d_%s', a, d{1});
                else
                    n = sprintf('KATP_ATP%d_%s_g%03d_Km%03d', a, d{1}, round(100 * g0Mult(gi)), round(100 * km));
                end
                y(gi) = getRow(n).dCaT_amp_pct;
            end
            styles = {'-o', '--s'};
            plot(ax, g0Mult, y, styles{1 + (km ~= 1)}, 'LineWidth', 1.5, 'Color', cols(j, :), 'MarkerFaceColor', cols(j, :));
            lab{end+1} = sprintf('ADP %.0f \\muM, K_m x%.3g', getRow(sprintf('ATP%d_%s', a, d{1})).ADP_uM, km); %#ok<SAGROW>
        end
    end
    yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
    set(ax, 'XTick', g0Mult); xlabel(ax, 'I_{KATP} g_0 multiplier'); ylabel(ax, '\Delta CaT amplitude vs 8 mM (%)');
    title(ax, sprintf('[ATP]_i = %d mM', a)); legend(ax, lab, 'Box', 'off', 'Location', 'southwest', 'FontSize', 8);
end
exportgraphics(f, fullfile(figDir, 'fig7_KATP_strength_band.png'), 'Resolution', 200);

% 8) Myofilament Ca2+ sensitivity (TnC ca50) and the free Ca transient
f = figure('Position', [100 100 560 400], 'Color', 'w'); ax = axes(f); hold(ax, 'on'); setStyle(ax);
xs = sort([1, ca50Mult]);
baseLab = mainLabels;
for b = 1:numel(baseSc)
    y = zeros(size(xs));
    for i = 1:numel(xs)
        if xs(i) == 1
            y(i) = getRow(baseSc{b}).dCaT_amp_pct;
        else
            y(i) = getRow(sprintf('ca50_%s_x%03d', baseSc{b}, round(100 * xs(i)))).dCaT_amp_pct;
        end
    end
    plot(ax, xs, y, '-o', 'LineWidth', 1.5, 'Color', cols(b, :), 'MarkerFaceColor', cols(b, :));
end
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
xlabel(ax, 'TnC ca50 multiplier (>1 = desensitised)'); ylabel(ax, '\Delta CaT amplitude vs 8 mM, ca50 x1 (%)');
legend(ax, baseLab, 'Box', 'off', 'Location', 'northwest');
title(ax, 'Myofilament Ca^{2+} sensitivity vs ATP effect');
exportgraphics(f, fullfile(figDir, 'fig8_ca50_sensitivity.png'), 'Resolution', 200);

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
