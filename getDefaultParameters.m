%     Cardiac model T-World
%     Copyright (C) 2026 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [parameters] = getDefaultParameters(parameters)
%GETDEFAULTPARAMETERS Populates the structure of parameters with defaults for those that are
%undefined.

% Cell type (endo/epi/mid)
if (~isfield(parameters,'cellType')) parameters.cellType = 0;end

% If the number of simulated beat is to be printed out.
if (~isfield(parameters,'verbose')) parameters.verbose = 0; end

% Extracellular concentrations
if (~isfield(parameters,'nao')) parameters.nao = 140; end
if (~isfield(parameters,'cao')) parameters.cao = 1.8; end
if (~isfield(parameters,'ko')) parameters.ko = 5; end


% Localization of ICaL and NCX: the fraction in junctional subspace
if (~isfield(parameters,'ICaL_fraction_junc')) parameters.ICaL_fraction_junc = 0.8; end
if (~isfield(parameters,'INaCa_fraction_junc')) parameters.INaCa_fraction_junc = 0.31; end

% Current multipliers
if (~isfield(parameters,'INa_Multiplier')) parameters.INa_Multiplier = 1; end
if (~isfield(parameters,'INaL_Multiplier')) parameters.INaL_Multiplier = 1; end
if (~isfield(parameters,'ICaL_Multiplier')) parameters.ICaL_Multiplier = 1; end
if (~isfield(parameters,'Itos_Multiplier')) parameters.Itos_Multiplier = 1; end
if (~isfield(parameters,'Itof_Multiplier')) parameters.Itof_Multiplier = 1; end
if (~isfield(parameters,'IKr_Multiplier')) parameters.IKr_Multiplier = 1; end
if (~isfield(parameters,'IKs_Multiplier')) parameters.IKs_Multiplier = 1; end
if (~isfield(parameters,'IK1_Multiplier')) parameters.IK1_Multiplier = 1; end
if (~isfield(parameters,'IKb_Multiplier')) parameters.IKb_Multiplier = 1; end
if (~isfield(parameters,'INaCa_Multiplier')) parameters.INaCa_Multiplier = 1; end
if (~isfield(parameters,'INaK_Multiplier')) parameters.INaK_Multiplier = 1; end
if (~isfield(parameters,'INab_Multiplier')) parameters.INab_Multiplier = 1; end
if (~isfield(parameters,'ICab_Multiplier')) parameters.ICab_Multiplier = 1; end
if (~isfield(parameters,'IpCa_Multiplier')) parameters.IpCa_Multiplier = 1; end
if (~isfield(parameters,'IClCa_Multiplier')) parameters.IClCa_Multiplier = 1; end
if (~isfield(parameters,'IClb_Multiplier')) parameters.IClb_Multiplier = 1; end
if (~isfield(parameters,'Jrel_Multiplier')) parameters.Jrel_Multiplier = 1; end
if (~isfield(parameters,'Jup_Multiplier')) parameters.Jup_Multiplier = 1; end

if (~isfield(parameters,'sex')) parameters.sex = 'U'; end % Unisex

% PKA signalling
if (~isfield(parameters,'PKA_P')) parameters.PKA_P = [0, 0, 0, 0, 0, 0, 0]; end 
if (~isfield(parameters,'Whole_cell_PP1')) parameters.Whole_cell_PP1 = 0.13698; end 
if (~isfield(parameters,'runSignallingPathway')) parameters.runSignallingPathway = 0; end 
if (~isfield(parameters,'ISO')) parameters.ISO = 0; end % Nonspecific

% Whole_cell_PP1
% PKA_P

% Intracellular nucleotides (mM), affecting SERCA, Na/K pump, sarcolemmal Ca pump
% (Cortassa et al. 2006) and I_KATP (Ferrero et al. 1996, human version of Mena Tobar
% et al. 2015). ATP_i and free ADP_i are held constant during a simulation (clamped).
% All effects are normalised to the reference (saturating) state ATP_ref/ADP_ref, so
% that ATP_i == ATP_ref and ADP_i == ADP_ref reproduce the unmodified model exactly.
% ADP_ref is the human creatine kinase equilibrium value at ATP_ref (getCKEquilibriumADP,
% ~94 uM). For reduced ATP, set ADP_i from measured or published free ADP: once
% phosphocreatine is exhausted the CK equilibrium no longer constrains ADP.
if (~isfield(parameters,'ATP_ref')) parameters.ATP_ref = 8.0; end
if (~isfield(parameters,'ADP_ref')) parameters.ADP_ref = getCKEquilibriumADP(parameters.ATP_ref); end
if (~isfield(parameters,'ATP_i')) parameters.ATP_i = parameters.ATP_ref; end
if (~isfield(parameters,'ADP_i')) parameters.ADP_i = parameters.ADP_ref; end
if (~isfield(parameters,'IKATP_Multiplier')) parameters.IKATP_Multiplier = 1; end % scales I_KATP maximal conductance g0
if (~isfield(parameters,'KATP_Km_Multiplier')) parameters.KATP_Km_Multiplier = 1; end % scales the K_ATP ATP half-inhibition Km (alpha)

% Myofilament (troponin C) Ca2+ sensitivity: multiplier of the Land model ca50 (>1
% desensitises, e.g. acidosis/Pi; <1 sensitises, e.g. MgADP). Changes TnC Ca2+
% buffering and hence the free Ca transient.
if (~isfield(parameters,'ca50_Multiplier')) parameters.ca50_Multiplier = 1; end

% An array of extra parameters (if user-defined, this can be a cell array
% as well), passing other parameters not defined otherwise in this script
if (~isfield(parameters,'extraParams')) parameters.extraParams = []; end

% There may be parameters defining clamp behaviour
if (~isfield(parameters,'vcParameters')) parameters.vcParameters = []; end
if (~isfield(parameters,'apClamp')) parameters.apClamp = []; end

% Stimulus parameters
if (~isfield(parameters,'stimAmp')) parameters.stimAmp = -53; end
if (~isfield(parameters,'stimDur')) parameters.stimDur = 1; end

% Model
if (~isfield(parameters,'model')) parameters.model = @model_TWorld; end


end

