% This script sweeps intracellular ATP ([ATP]_i, mM) and reports its effect on
% the Ca transient (amplitude, decay time constant) and the action potential.
% ATP and free ADP affect SERCA, the Na/K pump and the sarcolemmal Ca pump
% (Cortassa et al. 2006) and the ATP-sensitive K current I_KATP (Ferrero et al.
% 1996; human version by Mena Tobar et al. 2015). Free ADP is estimated from the
% creatine kinase equilibrium (getCKEquilibriumADP), given the residual
% phosphocreatine fraction and pH assumed at each ATP level.
% See analyses/lowATP for a fuller analysis and write-up.
%
% Note: lowering the Na/K pump rate makes [Na]_i drift over minutes, so values
% after a short prepacing are not a full steady state; increase 'beats' if needed.
clear

ATP_levels = [8.0, 4.0, 2.0]; % [mM]
fPCr_levels = [1.0, 0.10, 0.05]; % residual phosphocreatine (fraction of normoxic)
pH_levels = [7.1, 6.4, 6.2]; % intracellular pH
beats = 50; % beats simulated per ATP level; only the last is analysed
ignoreFirst = beats - 1;

X0 = getStartingState('TW_endo');
options = [];

results = table('Size', [numel(ATP_levels), 12], ...
    'VariableTypes', repmat({'double'}, 1, 12), ...
    'VariableNames', {'ATP_i', 'ADP_i_uM', 'fSERCA', 'fNaK', 'fKATP', 'CaT_amp_nM', 'CaT_dia_nM', 'CaT_tau_ms', 'APD90_ms', 'IKATP_peak', 'CaSR_mM', 'Nai_mM'});

figure(1); clf;
for iATP = 1:numel(ATP_levels)
    param = [];
    param.model = @model_TWorld;
    param.bcl = 1000;
    param.ATP_i = ATP_levels(iATP);
    param.ADP_i = getCKEquilibriumADP(ATP_levels(iATP), fPCr_levels(iATP), pH_levels(iATP));

    [time, X] = modelRunner(X0, options, param, beats, ignoreFirst);
    currents = getCurrentsStructure(time, X, param, 0);

    t = currents.time;
    cai = currents.Ca_i * 1e6; % nM
    [caPeak, iPeak] = max(cai);
    caDia = min(cai);

    results.ATP_i(iATP) = ATP_levels(iATP);
    results.ADP_i_uM(iATP) = 1000 * param.ADP_i;
    results.fSERCA(iATP) = currents.fATP_SERCA(1);
    results.fNaK(iATP) = currents.fATP_NaK(1);
    results.fKATP(iATP) = currents.fATP_KATP(1);
    results.CaT_amp_nM(iATP) = caPeak - caDia;
    results.CaT_dia_nM(iATP) = caDia;
    results.CaT_tau_ms(iATP) = getDecayTau(t, cai, iPeak);
    results.APD90_ms(iATP) = getAPD(t, currents.V, 0.9);
    results.IKATP_peak(iATP) = max(currents.IKATP);
    results.CaSR_mM(iATP) = max(currents.Ca_SR);
    results.Nai_mM(iATP) = currents.Na_i(end);

    subplot(1, 2, 1); hold on; plot(t, currents.V);
    subplot(1, 2, 2); hold on; plot(t, cai);
end
subplot(1, 2, 1); xlabel('Time (ms)'); ylabel('Membrane potential (mV)');
legend(arrayfun(@(a) sprintf('ATP_i = %g mM', a), ATP_levels, 'UniformOutput', false));
subplot(1, 2, 2); xlabel('Time (ms)'); ylabel('[Ca]_i (nM)');

disp(results);
ref = results(1, :);
for iATP = 2:height(results)
    fprintf('ATP_i = %4.2f mM vs %4.2f mM: CaT amplitude %+6.1f %%, decay tau %+6.1f %%, APD90 %+6.1f %%\n', ...
        results.ATP_i(iATP), ref.ATP_i, ...
        100 * (results.CaT_amp_nM(iATP) / ref.CaT_amp_nM - 1), ...
        100 * (results.CaT_tau_ms(iATP) / ref.CaT_tau_ms - 1), ...
        100 * (results.APD90_ms(iATP) / ref.APD90_ms - 1));
end

%% Local functions
function tau = getDecayTau(t, cai, iPeak)
% Single-exponential fit, cai = A*exp(-(t-t0)/tau) + C, to the decay phase
% starting when 10% of the transient amplitude has decayed (excluding the peak).
caDia = min(cai(iPeak:end));
amp = cai(iPeak) - caDia;
iStart = iPeak - 1 + find(cai(iPeak:end) <= caDia + 0.9 * amp, 1);
tFit = t(iStart:end) - t(iStart);
cFit = cai(iStart:end);
sse = @(p) sum((p(1) * exp(-tFit / p(2)) + p(3) - cFit).^2);
p = fminsearch(sse, [cFit(1) - cFit(end), 100, cFit(end)], optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4));
tau = p(2);
end

function apd = getAPD(t, v, fraction)
% APD at the given repolarisation fraction (e.g. 0.9 for APD90).
[vPeak, iPeak] = max(v);
vRest = v(1);
vThresh = vPeak - fraction * (vPeak - vRest);
iRepol = iPeak - 1 + find(v(iPeak:end) <= vThresh, 1);
iUp = find(v >= vThresh, 1);
apd = interp1(v(iRepol-1:iRepol), t(iRepol-1:iRepol), vThresh) - interp1(v(iUp-1:iUp), t(iUp-1:iUp), vThresh);
end
