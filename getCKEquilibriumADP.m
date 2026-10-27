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

function [ADP, state] = getCKEquilibriumADP(ATP, fPCr, pH)
%GETCKEQUILIBRIUMADP Free cytosolic [ADP] (mM) from the creatine kinase equilibrium.
%
% PCr + ADP + H+ <-> Cr + ATP is assumed to be at equilibrium (CK activity in the
% heart far exceeds ATP turnover), so
%   [ADP] = [ATP][Cr] / (K'CK(pH) [PCr]),   K'CK(pH) = 166 * 10^(-0.87 (pH - 7)),
% K'CK = 166 at pH 7.0 (Lawson & Veech 1979, 1.66e9 M^-1 at 38 degC, free Mg 1 mM),
% with the fractional proton stoichiometry of 0.87 (as used by Stepanov et al. 1997).
%
% INPUTS:
% ATP - cytosolic [ATP] (mM)
% fPCr - [PCr] as a fraction of the normoxic [PCr] (default 1)
% pH - intracellular pH (default 7.1)
%
% The total creatine pool is conserved; normoxic pools are referenced to [ATP] = 8 mM:
% PCr/ATP = 1.6 (human cardiac 31P-MRS) and PCr/total creatine = 0.6.
%
% OUTPUTS:
% ADP - free [ADP] (mM)
% state - structure with ATP, ADP, PCr, Cr, Cr_tot, pH and KCK

if nargin < 2 || isempty(fPCr), fPCr = 1; end
if nargin < 3 || isempty(pH), pH = 7.1; end

ATP_normoxic = 8;                 % [mM]
PCr_normoxic = 1.6 * ATP_normoxic; % [mM]
Cr_tot = PCr_normoxic / 0.6;      % [mM]

PCr = fPCr * PCr_normoxic;
Cr = Cr_tot - PCr;
KCK = 166 * 10^(-0.87 * (pH - 7));
ADP = ATP * Cr / (KCK * PCr);

state = struct('ATP', ATP, 'ADP', ADP, 'PCr', PCr, 'Cr', Cr, 'Cr_tot', Cr_tot, 'pH', pH, 'KCK', KCK);
end
