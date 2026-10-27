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
%   [ADP] = [ATP][Cr] / (K_CK [H+] [PCr]),   K_CK = 1.66e9 M^-1
% (Lawson & Veech 1979, 38 degC, free Mg 1 mM; K' = 166 at pH 7.0), as used for the
% human heart by Weiss et al. 2005 (PNAS 102:808).
%
% INPUTS:
% ATP - cytosolic [ATP] (mM)
% fPCr - [PCr] as a fraction of the normoxic [PCr] (default 1)
% pH - intracellular pH (default 7.05)
%
% The total creatine pool is conserved; normoxic pools are referenced to [ATP] = 8 mM
% using healthy human LV 31P/1H-MRS (Weiss et al. 2005): PCr/ATP = 10.1/5.7 and total
% creatine/ATP = 27.6/5.7 (umol/g wet weight), giving free ADP ~ 94 uM at 8 mM ATP.
% The equilibrium is ill-conditioned when PCr is depleted (ischaemia), where measured
% or published free ADP values should be used instead.
%
% OUTPUTS:
% ADP - free [ADP] (mM)
% state - structure with ATP, ADP, PCr, Cr, Cr_tot, pH and KCK

if nargin < 2 || isempty(fPCr), fPCr = 1; end
if nargin < 3 || isempty(pH), pH = 7.05; end

ATP_normoxic = 8;                         % [mM]
PCr_normoxic = 10.1 / 5.7 * ATP_normoxic; % [mM]
Cr_tot = 27.6 / 5.7 * ATP_normoxic;       % [mM]

PCr = fPCr * PCr_normoxic;
Cr = Cr_tot - PCr;
KCK = 166 * 10^(-(pH - 7));
ADP = ATP * Cr / (KCK * PCr);

state = struct('ATP', ATP, 'ADP', ADP, 'PCr', PCr, 'Cr', Cr, 'Cr_tot', Cr_tot, 'pH', pH, 'KCK', KCK);
end
