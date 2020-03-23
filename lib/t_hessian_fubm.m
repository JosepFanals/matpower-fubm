function t_hessian_fubm(quiet)
%T_HESSIAN_FUBM  Numerical tests of 2nd derivative code.

%   This code compares the results from the obtained derivatives against
%   the aproximated derivatives using the finite differences method.

%   FINITE DIFFERENCES METHOD
%   This method calculates the derivatives with an aproximation as:
%   f' (x) ~~ ( f(x+h) - f(x) ) / h 
%   f''(x) ~~ ( f'(x+h) - f'(x) ) / h 

%   ABRAHAM ALVAREZ BUSTOS
%   This code is based and created for MATPOWER
%   This is part of the Flexible Universal Branch Model (FUBM) for Matpower
%   For more info about the model, email: 
%   snoop_and@hotmail.com, abraham.alvarez-bustos@durham.ac.uk 
%
%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(96, quiet); %AAB-initializes the global test counters (Number of Total Tests)
%casefile = 'case30';
%casefile = 'fubm_caseHVDC_qt';
%casefile = 'fubm_caseHVDC_vt';
%casefile = 'fubm_case_57_14_2MTDC_ctrls';
%casefile = 'fubm_case_30_2MTDC_ctrls_vt1_pf';
casefile = 'fubm_case_30_2MTDC_ctrls_vt2_pf';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX, VF_SET, VT_SET,TAP_MAX, ...
    TAP_MIN, CONV, BEQ, K2, BEQ_MIN, BEQ_MAX, SH_MIN, SH_MAX, GSW, ...
    ALPH1, ALPH2, ALPH3] = idx_brch;%<<AAB-extra fields for FUBM- Original: idx_brch
vcart=0;
%% run powerflow to get solved case
mpopt = mpoption('verbose', 0, 'out.all', 0);
[baseMVA, bus, gen, branch, success, et] = runpf(casefile, mpopt);

%% switch to internal bus numbering and build admittance matrices
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.
Vm = bus(:, VM);
Va = bus(:, VA) * pi/180;
V = Vm .* exp(1j * Va);
Vr = real(V);
Vi = imag(V);
f = branch(:, F_BUS);       %% list of "from" buses
t = branch(:, T_BUS);       %% list of "to" buses
nl = length(f);
nb = length(V);
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
pert = 1e-8;                                    %% perturbation factor (h) for the Finite Differences Method

%% identifier of AC/DC grids
iBeqz = find (branch(:,CONV)==1 & branch(:, BR_STATUS)==1); %AAB- Find branch locations of VSC, If the grid has them it's an AC/DC grid
nBeqz = length(iBeqz); %AAB- Number of VSC with active Zero Constraint control
%%identifier of elements with Vf controlled by Beq
iBeqv = find (branch(:,CONV)==2 & branch(:, BR_STATUS)==1 & branch(:, VF_SET)~=0); %AAB- Find branch locations of VSC
if nBeqz
    nBeqv = length(iBeqv); %AAB- Number of VSC with Vf controlled by Beq
else
    nBeqv = 0; %AAB- Vdc control with Beq requires an AC/DC grid.
end
iVscL = find (branch(:,CONV)~=0 & branch(:, BR_STATUS)==1 & (branch(:, ALPH1)~=0 | branch(:, ALPH2)~=0 | branch(:, ALPH3)~=0) ); %AAB- Find VSC with active PWM Losses Calculation [nVscL,1]
if nBeqz
    nVscL = length(iVscL); %AAB- Number of VSC with power losses
else
    nVscL = 0; %AAB- Number of VSC with power losses
end
%% Identify if grid has controls
iPfsh = find (branch(:,PF)~=0 & branch(:, BR_STATUS)==1 & (branch(:, SH_MIN)~=-360 | branch(:, SH_MAX)~=360)); %AAB- Find branch locations with Pf controlled by Theta_shift [nPfsh,1]
nPfsh = length(iPfsh); %AAB- Number of elements with active Pf controlled by Theta_shift
iQtma = find (branch(:,QT)~=0 &branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:,VT_SET)==0 ); %AAB- Find branch locations with Qt controlled by ma/tap [nQtma,1]
nQtma = length(iQtma); %AAB- Number of elements with active Qt controlled by ma/tap
iVtma = find (branch(:, BR_STATUS)==1 & (branch(:, TAP_MIN)~= branch(:, TAP_MAX)) & branch(:, VT_SET)~=0 ); %AAB- Find branch locations with Vt controlled by ma/tap [nVtma,1]
nVtma = length(iVtma); %AAB- Number of elements with active Vt controlled by ma/tap

%% -----  run tests for polar coordinate  -----
    %%-----  create perturbed voltages  -----
    %% polar coordinate voltages (V1=Va, V2=Vm)
        coord = 'polar';
        vv = {'aa', 'av', 'va', 'vv'};
        V1p = (Vm*ones(1,nb)) .* (exp(1j * (Va*ones(1,nb) + pert*eye(nb,nb))));
        V2p = (Vm*ones(1,nb) + pert*eye(nb,nb)) .* (exp(1j * Va) * ones(1,nb));

    %% -----  check d2Sbus_dV2 code  -----
    t = ' - d2Sbus_dV2 (complex power injections)';
    lam = 10 * rand(nb, 1);
    %%sparse matrices partial derivatives
    [H11, H12, H21, H22] = d2Sbus_dV2(Ybus, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    num_H11 = zeros(nb, nb);
    num_H12 = zeros(nb, nb);
    num_H21 = zeros(nb, nb);
    num_H22 = zeros(nb, nb);
    [dSbus_dV1, dSbus_dV2] = dSbus_dV(Ybus, V, vcart); %dSbus_dVa
    %VaVa
    for i = 1:nb
        V1p = V;
        V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
        [dSbus_dV1_1p, dSbus_dV2_1p] = dSbus_dV(Ybus, V1p, vcart); %dSbus_dVaPertVa
        num_H11(:, i) = (dSbus_dV1_1p - dSbus_dV1).' * lam / pert; %VaVa (dSbus_dVaPertVa - dSbus_dVa)
    end
    %VaVm
    for i = 1:nb
        V2p = V;
        V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSbus_dV1_2p, dSbus_dV2_2p] = dSbus_dV(Ybus, V2p, vcart); %dSbus_dVaPertVm
        num_H12(:, i) = (dSbus_dV1_2p - dSbus_dV1).' * lam / pert; %VaVm (dSbus_dVaPertVm - dSbus_dVa)
    end
    %VmVa
    for i = 1:nb
        V1p = V;
        V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
        [dSbus_dV1_1p, dSbus_dV2_1p] = dSbus_dV(Ybus, V1p, vcart); %dSbus_dVmPertVa
        num_H21(:, i) = (dSbus_dV2_1p - dSbus_dV2).' * lam / pert; %VmVa (dSbus_dVmPertVa - dSbus_dVm)
    end
    %VmVm
    for i = 1:nb
        V2p = V;
        V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSbus_dV1_2p, dSbus_dV2_2p] = dSbus_dV(Ybus, V2p, vcart); %dSbus_dVmPertVm
        num_H22(:, i) = (dSbus_dV2_2p - dSbus_dV2).' * lam / pert; %VmVm (dSbus_dVmPertVm - dSbus_dVm)
    end
    
    t_is(full(H11), num_H11, 4, sprintf('%s - H%s%s', coord, vv{1}, t));
    t_is(full(H12), num_H12, 4, sprintf('%s - H%s%s', coord, vv{2}, t));
    
    t_is(full(H21), num_H21, 4, sprintf('%s - H%s%s', coord, vv{3}, t));
    t_is(full(H22), num_H22, 4, sprintf('%s - H%s%s', coord, vv{4}, t));


    %% -----  check d2Sbus_dxBeqz2 code  -----
    t = ' - d2Sbus_dxBeqz2 (Beqz complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G15, G25, G51, G52, G55] = d2Sbus_dxBeqz2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G15, num_G25, num_G51, num_G52, num_G55] = d2Sbus_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G15), num_G15, 4, sprintf('%s - HVaBeqz%s', coord, t));
    t_is(full(G25), num_G25, 4, sprintf('%s - HVmBeqz%s', coord, t));
    
    t_is(full(G51), num_G51, 4, sprintf('%s - HBeqzVa%s', coord, t));
    t_is(full(G52), num_G52, 4, sprintf('%s - HBeqzVm%s', coord, t));
    
    t_is(full(G55), num_G55, 4, sprintf('%s - HBeqz2 %s', coord, t));

    %% -----  check d2Sbus_dxBeqv2 code  -----
    t = ' - d2Sbus_dxBeqv2 (Beqv complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G16, G26, G56, G61, G62, G65, G66] = d2Sbus_dxBeqv2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G16, num_G26, num_G56, num_G61, num_G62, num_G65, num_G66] = d2Sbus_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G16), num_G16, 4, sprintf('%s - HVaBeqv%s'  , coord, t));
    t_is(full(G26), num_G26, 4, sprintf('%s - HVmBeqv%s'  , coord, t));
    t_is(full(G56), num_G56, 4, sprintf('%s - HBeqzBeqv%s', coord, t));
    
    t_is(full(G61), num_G61, 4, sprintf('%s - HBeqvVa%s'  , coord, t));
    t_is(full(G62), num_G62, 4, sprintf('%s - HBeqvVm%s'  , coord, t));
    t_is(full(G65), num_G65, 4, sprintf('%s - HBeqvBeqz%s', coord, t));
    
    t_is(full(G66), num_G66, 4, sprintf('%s - HBeqv2 %s'  , coord, t));

    %% -----  check d2Sbus_dxsh2 code  -----
    t = ' - d2Sbus_dxsh2 (Pf_sh complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G17, G27, G57, G67, G71, G72, G75, G76, G77] = d2Sbus_dxsh2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G17, num_G27, num_G57, num_G67, num_G71, num_G72, num_G75, num_G76, num_G77] = d2Sbus_dxsh2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G17), num_G17, 4, sprintf('%s - HVaPfsh%s'  , coord, t));
    t_is(full(G27), num_G27, 4, sprintf('%s - HVmPfsh%s'  , coord, t));
    t_is(full(G57), num_G57, 4, sprintf('%s - HBeqzPfsh%s', coord, t));
    t_is(full(G67), num_G67, 4, sprintf('%s - HBeqvPfsh%s', coord, t));
    
    t_is(full(G71), num_G71, 4, sprintf('%s - HPfshVa%s'  , coord, t));
    t_is(full(G72), num_G72, 4, sprintf('%s - HPfshVm%s'  , coord, t));
    t_is(full(G75), num_G75, 4, sprintf('%s - HPfshBeqz%s', coord, t));
    t_is(full(G76), num_G76, 4, sprintf('%s - HPfshBeqv%s', coord, t));
    
    t_is(full(G77), num_G77, 4, sprintf('%s - HPfsh2 %s'  , coord, t));

    %% -----  check d2Sbus_dxqtma2 code  -----
    t = ' - d2Sbus_dxqtma2 (Qt_ma complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G18, G28, G58, G68, G78, G81, G82, G85, G86, G87, G88] = d2Sbus_dxqtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G18, num_G28, num_G58, num_G68, num_G78, num_G81, num_G82, num_G85, num_G86, num_G87, num_G88] = d2Sbus_dxqtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G18), num_G18, 4, sprintf('%s - HVaQtma%s'  , coord, t));
    t_is(full(G28), num_G28, 4, sprintf('%s - HVmQtma%s'  , coord, t));
    t_is(full(G58), num_G58, 4, sprintf('%s - HBeqzQtma%s', coord, t));
    t_is(full(G68), num_G68, 4, sprintf('%s - HBeqvQtma%s', coord, t));
    t_is(full(G78), num_G78, 4, sprintf('%s - HPfshQtma%s', coord, t));
    
    t_is(full(G81), num_G81, 4, sprintf('%s - HQtmaVa%s'  , coord, t));
    t_is(full(G82), num_G82, 4, sprintf('%s - HQtmaVm%s'  , coord, t));
    t_is(full(G85), num_G85, 4, sprintf('%s - HQtmaBeqz%s', coord, t));
    t_is(full(G86), num_G86, 4, sprintf('%s - HQtmaBeqv%s', coord, t));
    t_is(full(G87), num_G87, 4, sprintf('%s - HQtmaPfsh%s', coord, t));
    
    t_is(full(G88), num_G88, 4, sprintf('%s - HQtma2 %s'  , coord, t));
    
    %% -----  check d2Sbus_dxvtma2 code  -----
    t = ' - d2Sbus_dxvtma2 (Vt_ma complex power injections)';
    lam = 10 * rand(nb   , 1    );
    %%sparse matrices partial derivatives
    [G19, G29, G59, G69, G79, G89, G91, G92, G95, G96, G97, G98, G99] = d2Sbus_dxvtma2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_G19, num_G29, num_G59, num_G69, num_G79, num_G89, num_G91, num_G92, num_G95, num_G96, num_G97, num_G98, num_G99] = d2Sbus_dxvtma2Pert(baseMVA, bus, branch, V, lam, pert, vcart);  

    t_is(full(G19), num_G19, 4, sprintf('%s - HVaVtma%s'  , coord, t));
    t_is(full(G29), num_G29, 4, sprintf('%s - HVmVtma%s'  , coord, t));
    t_is(full(G59), num_G59, 4, sprintf('%s - HBeqzVtma%s', coord, t));
    t_is(full(G69), num_G69, 4, sprintf('%s - HBeqvVtma%s', coord, t));
    t_is(full(G79), num_G79, 4, sprintf('%s - HPfshVtma%s', coord, t));
    t_is(full(G89), num_G89, 4, sprintf('%s - HQtmaVtma%s', coord, t));
    
    t_is(full(G91), num_G91, 4, sprintf('%s - HVtmaVa%s'  , coord, t));
    t_is(full(G92), num_G92, 4, sprintf('%s - HVtmaVm%s'  , coord, t));
    t_is(full(G95), num_G95, 4, sprintf('%s - HVtmaBeqz%s', coord, t));
    t_is(full(G96), num_G96, 4, sprintf('%s - HVtmaBeqv%s', coord, t));
    t_is(full(G97), num_G97, 4, sprintf('%s - HVtmaPfsh%s', coord, t));
    t_is(full(G98), num_G98, 4, sprintf('%s - HVtmaQtma%s', coord, t));
    
    t_is(full(G99), num_G99, 4, sprintf('%s - HVtma2 %s'  , coord, t));

    %% -----  check d2Sbr_dV2 code  -----
    t = ' - d2Sbr_dV2 (complex power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [Gf11, Gf12, Gf21, Gf22] = d2Sbr_dV2(Cf, Yf, V, lam, vcart);
    [Gt11, Gt12, Gt21, Gt22] = d2Sbr_dV2(Ct, Yt, V, lam, vcart);
    for i = 1:nb
        V1p = V;
        V2p = V;
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        num_Gf11(:, i) = (dSf_dV1_1p - dSf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dSf_dV2_1p - dSf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dSt_dV1_1p - dSt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dSt_dV2_1p - dSt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        num_Gf12(:, i) = (dSf_dV1_2p - dSf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dSf_dV2_2p - dSf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dSt_dV1_2p - dSt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dSt_dV2_2p - dSt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 4, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 4, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 4, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 4, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 4, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 4, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 4, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 4, sprintf('%s - Gt%s%s', coord, vv{4}, t));

    %% -----  check d2Sbr_dxBeqz2 code  -----
    t = ' - d2Sbr_dxBeqz2 (Beqz complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf13, Hf23, Hf31, Hf32, Hf33] = d2Sf_dxBeqz2(branch, V, lam, vcart);
    [Ht13, Ht23, Ht31, Ht32, Ht33] = d2St_dxBeqz2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf13, num_Hf23, num_Hf31, num_Hf32, num_Hf33,...
     num_Ht13, num_Ht23, num_Ht31, num_Ht32, num_Ht33] = d2Sbr_dxBeqz2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf13), num_Hf13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Hf23), num_Hf23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    
    t_is(full(Hf31), num_Hf31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Hf32), num_Hf32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    
    t_is(full(Hf33), num_Hf33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht13), num_Ht13, 4, sprintf('%s - HVaBeqz%s%s', coord, t, br));
    t_is(full(Ht23), num_Ht23, 4, sprintf('%s - HVmBeqz%s%s', coord, t, br));
    
    t_is(full(Ht31), num_Ht31, 4, sprintf('%s - HBeqzVa%s%s', coord, t, br));
    t_is(full(Ht32), num_Ht32, 4, sprintf('%s - HBeqzVm%s%s', coord, t, br));
    
    t_is(full(Ht33), num_Ht33, 4, sprintf('%s - HBeqz2 %s%s', coord, t, br));

    %% -----  check d2Sbr_dxBeqv2 code  -----
    t = ' - d2Sbr_dxBeqv2 (Beqv complex power flows)';
    lam = 10 * rand(nl, 1);
    %%sparse matrices partial derivatives
    [Hf14, Hf24, Hf34, Hf41, Hf42, Hf43, Hf44] = d2Sf_dxBeqv2(branch, V, lam, vcart);
    [Ht14, Ht24, Ht34, Ht41, Ht42, Ht43, Ht44] = d2St_dxBeqv2(branch, V, lam, vcart);
    
    %%compute numerically to compare (Finite Differences Method)
    [num_Hf14, num_Hf24, num_Hf34, num_Hf41, num_Hf42, num_Hf43, num_Hf44,...
     num_Ht14, num_Ht24, num_Ht34, num_Ht41, num_Ht42, num_Ht43, num_Ht44] = d2Sbr_dxBeqv2Pert(baseMVA, bus, branch, V, lam, pert, vcart);
    
    br = ' - "from" side';
    t_is(full(Hf14), num_Hf14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Hf24), num_Hf24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Hf34), num_Hf34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    
    t_is(full(Hf41), num_Hf41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Hf42), num_Hf42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Hf43), num_Hf43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    
    t_is(full(Hf44), num_Hf44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    br = ' - " to " side';
    t_is(full(Ht14), num_Ht14, 4, sprintf('%s - HVaBeqv%s%s', coord, t, br));
    t_is(full(Ht24), num_Ht24, 4, sprintf('%s - HVmBeqv%s%s', coord, t, br));
    t_is(full(Ht34), num_Ht34, 4, sprintf('%s - HBeqzBeqv%s%s', coord, t, br));
    
    t_is(full(Ht41), num_Ht41, 4, sprintf('%s - HBeqvVa%s%s', coord, t, br));
    t_is(full(Ht42), num_Ht42, 4, sprintf('%s - HBeqvVm%s%s', coord, t, br));
    t_is(full(Ht43), num_Ht43, 4, sprintf('%s - HBeqvBeqz%s%s', coord, t, br));
    
    t_is(full(Ht44), num_Ht44, 4, sprintf('%s - HBeqv2 %s%s', coord, t, br));
    
    
    %% -----  check d2Abr_dV2 code  -----
    t = ' - d2Abr_dV2 (squared apparent power flows)';
    lam = 10 * rand(nl, 1);
    % lam = [1; zeros(nl-1, 1)];
    num_Gf11 = zeros(nb, nb);
    num_Gf12 = zeros(nb, nb);
    num_Gf21 = zeros(nb, nb);
    num_Gf22 = zeros(nb, nb);
    num_Gt11 = zeros(nb, nb);
    num_Gt12 = zeros(nb, nb);
    num_Gt21 = zeros(nb, nb);
    num_Gt22 = zeros(nb, nb);
    d2Sf_dV2 = @(V, mu)d2Sbr_dV2(Cf, Yf, V, mu, vcart);
    d2St_dV2 = @(V, mu)d2Sbr_dV2(Ct, Yt, V, mu, vcart);
    [dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St] = dSbr_dV(branch, Yf, Yt, V, vcart);
    [dAf_dV1, dAf_dV2, dAt_dV1, dAt_dV2] = ...
                            dAbr_dV(dSf_dV1, dSf_dV2, dSt_dV1, dSt_dV2, Sf, St);
    [Gf11, Gf12, Gf21, Gf22] = d2Abr_dV2(d2Sf_dV2, dSf_dV1, dSf_dV2, Sf, V, lam);
    [Gt11, Gt12, Gt21, Gt22] = d2Abr_dV2(d2St_dV2, dSt_dV1, dSt_dV2, St, V, lam);
    for i = 1:nb
        V1p = V;
        V2p = V;
            V1p(i) = Vm(i) * exp(1j * (Va(i) + pert));  %% perturb Va
            V2p(i) = (Vm(i) + pert) * exp(1j * Va(i));  %% perturb Vm
        [dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p] = ...
            dSbr_dV(branch, Yf, Yt, V1p, vcart);
        [dAf_dV1_1p, dAf_dV2_1p, dAt_dV1_1p, dAt_dV2_1p] = ...
            dAbr_dV(dSf_dV1_1p, dSf_dV2_1p, dSt_dV1_1p, dSt_dV2_1p, Sf_1p, St_1p);
        num_Gf11(:, i) = (dAf_dV1_1p - dAf_dV1).' * lam / pert;
        num_Gf21(:, i) = (dAf_dV2_1p - dAf_dV2).' * lam / pert;
        num_Gt11(:, i) = (dAt_dV1_1p - dAt_dV1).' * lam / pert;
        num_Gt21(:, i) = (dAt_dV2_1p - dAt_dV2).' * lam / pert;

        [dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p] = ...
            dSbr_dV(branch, Yf, Yt, V2p, vcart);
        [dAf_dV1_2p, dAf_dV2_2p, dAt_dV1_2p, dAt_dV2_2p] = ...
            dAbr_dV(dSf_dV1_2p, dSf_dV2_2p, dSt_dV1_2p, dSt_dV2_2p, Sf_2p, St_2p);
        num_Gf12(:, i) = (dAf_dV1_2p - dAf_dV1).' * lam / pert;
        num_Gf22(:, i) = (dAf_dV2_2p - dAf_dV2).' * lam / pert;
        num_Gt12(:, i) = (dAt_dV1_2p - dAt_dV1).' * lam / pert;
        num_Gt22(:, i) = (dAt_dV2_2p - dAt_dV2).' * lam / pert;
    end

    t_is(full(Gf11), num_Gf11, 2.5, sprintf('%s - Gf%s%s', coord, vv{1}, t));
    t_is(full(Gf12), num_Gf12, 2.5, sprintf('%s - Gf%s%s', coord, vv{2}, t));
    t_is(full(Gf21), num_Gf21, 2.5, sprintf('%s - Gf%s%s', coord, vv{3}, t));
    t_is(full(Gf22), num_Gf22, 2.5, sprintf('%s - Gf%s%s', coord, vv{4}, t));

    t_is(full(Gt11), num_Gt11, 2.5, sprintf('%s - Gt%s%s', coord, vv{1}, t));
    t_is(full(Gt12), num_Gt12, 2.5, sprintf('%s - Gt%s%s', coord, vv{2}, t));
    t_is(full(Gt21), num_Gt21, 2.5, sprintf('%s - Gt%s%s', coord, vv{3}, t));
    t_is(full(Gt22), num_Gt22, 2.5, sprintf('%s - Gt%s%s', coord, vv{4}, t));

t_end;
