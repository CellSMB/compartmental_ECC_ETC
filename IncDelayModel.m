
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = IncDelayModel(params,inits)
    fluxvar   =params(1);
    t_max     =params(2);
    K_c       =params(3);
    K_h       =params(4);
    K_t       =params(5);
    global var1;
    var1=fluxvar;
    global var3;
    var3 = t_max;
    global var4;
    var4 = K_c;
    global var5;
    var5 = K_h;
    global var6;
    var6 = K_t;
    global initStates;
    initStates=inits;
    % [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel([0.25,1,1,3,10,1.5]);
% [VOI2, STATES2, ALGEBRAIC2, CONSTANTS2]=IncDelayModel([1,1,1,1,1,1]);
% [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel([1,1,5000,1,1,1]);
%[VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel([0.25,1,10000,3,10,1.5]);
% [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel([1,1,1,2,3]);
 % [VOI, STATES, ALGEBRAIC, CONSTANTS]=IncDelayModel([1,1,10000,2,3,1]);
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
%    exp = max(1-ALGEBRAIC(:,59).*ALGEBRAIC(:,60))
%    max(ALGEBRAIC(:,66))
%    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
% trapz(VOI(time),ALGEBRAIC(time,66))
%    plotSimulation(VOI, STATES, ALGEBRAIC, CONSTANTS)
end

function plotSimulation(VOI, STATES, ALGEBRAIC, CONSTANTS)
%% Plot h's
global var3;
if 1
    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    figure('pos',[600,0,600,600])
    numa = 1;
    numd =3;
    subplot(numd,numa,1)
    hold on
%     plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
    legend({'Ca_i','Ca_i+IP_3'})
%     title(var3);
    hold on
    subplot(numd,numa,2)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
    legend('IP_3R')
    subplot(numd,numa,3)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,65))
    legend('P_0')
end
if 0
    time=VOI>(VOI(end)-CONSTANTS(:,7)-500);
    figure('pos',[0,0,600,600])
    numa = 1;
    numd =3;
    subplot(numd,numa,1)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
    plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
    legend({'Ca_i+IP_3','Ca_i'})
    title(var3);
    hold on
    subplot(numd,numa,2)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,7))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,60))
    legend({'h','h_{alpha}'})
    subplot(numd,numa,3)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
end
if 0
   figure('pos',[600,0,600,600])
    numa = 1;
    numd =4;
    subplot(numd,numa,1)
    hold on
    plot(VOI(1:(CONSTANTS(:,7)))+500,STATES(1:(CONSTANTS(:,7)),1))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),STATES(time,1))
    legend({'Ca_i','Ca_i+IP_3'})
    title(var3);
    hold on
    subplot(numd,numa,2)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,66))
    legend('IP_3R')
    subplot(numd,numa,3)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,63))
    legend('alpha')
    subplot(numd,numa,4)
    hold on
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)-500),ALGEBRAIC(time,64))
    legend('beta')
end
if 1 
    figure('pos',[1200,0,600,400])
    ca = 0:1e-5:1e-3;
    plot(ca,power(ca, 4.00000)./(power(CONSTANTS(:,50), 4.00000)+power(ca, 4.00000)))
    m=power(ca, 4.00000)./(power(CONSTANTS(:,50), 4.00000)+power(ca, 4.00000));
    h_alpha=power(CONSTANTS(:,51), 4.00000)./(power(CONSTANTS(:,51), 4.00000)+power(ca, 4.00000));
    hold on
    plot(ca,h_alpha)
    plot(ca,m.*h_alpha)
end
if 0
    time=VOI>(VOI(end)-CONSTANTS(:,7));
    figure
    hold on
    plot(VOI(1:CONSTANTS(:,7)),STATES(1:CONSTANTS(:,7),1))
    plot(VOI(time)-(VOI(end)-CONSTANTS(:,7)),STATES(time,1))
end
%% Plot calcium fluxes
    if 0
        figure
        numa = 2;
        numd =3;
        start = 1;%floor(1.8*CONSTANTS(:,7));
        subplot(numd,numa,1)
        hold on
        plot(VOI(start:end),STATES(start:end,1))
        title('Ca_{cyto}')
        subplot(numd,numa,2)
        hold on
        plot(VOI(start:end),STATES(start:end,2))
        title('Ca_{nSR}')
        
        subplot(numd,numa,numa+1)
        hold on
        plot(VOI(start:end),ALGEBRAIC(start:end,35))
        title('I_{RyR}')
        subplot(numd,numa,numa+2)
        hold on
        plot(VOI(start:end),ALGEBRAIC(start:end,44))
        title('I_{LCC}')
        
        subplot(numd,numa,2*numa+1)
        hold on
        plot(VOI(start:end),ALGEBRAIC(start:end,66))
        title('I_{IP3R}')
        subplot(numd,numa,2*numa+2)
        hold on
        plot(VOI(start:end),ALGEBRAIC(start:end,48))
        title('I_{SERCA}')
        
%         subplot(numd,numa,3*numa+1)
%         hold on
%         
%         subplot(numd,numa,3*numa+2)
%         hold on
    end
%% Plot just Ca_i, Ca_jSR, Ca_nSR
    if 0
        figure
        start=1000;
        subplot(2,1,1)
        hold on
        plot(VOI(start:end),STATES(start:end,1))
        title('Ca_{cyto}')
        subplot(2,1,2)
        hold on
        plot(VOI(start:end),STATES(start:end,2))
        title('Ca_{SR}')
    end
    
%% Plot all fluxexs in and out of cytosol
if 0
    figure
    numa = 3;
    numd = 5;
    start = 1;%floor(1.8*CONSTANTS(:,7));
    subplot(numd,numa,1)
    hold on
    plot(VOI(start:end),STATES(start:end,1))
    title('Ca_{cyto}')
    subplot(numd,numa,2)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,58))
    title('I_{other}')
    subplot(numd,numa,3)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,54))
    title('I_{SR}')

    subplot(numd,numa,numa+1)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,35))
    title('I_{RyR}')
    subplot(numd,numa,numa+2)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,44))
    title('I_{LCC}')
    subplot(numd,numa,numa+3)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,46))
    title('I_{NaCa}')

    subplot(numd,numa,2*numa+1)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,66))
    title('I_{IP3R}')
    subplot(numd,numa,2*numa+2)
    hold on
    plot(VOI(start:end),-ALGEBRAIC(start:end,48))
    title('I_{SERCA}')
    subplot(numd,numa,2*numa+3)
    hold on
    plot(VOI(start:end),-ALGEBRAIC(start:end,50))
    title('I_{pCa}')

    subplot(numd,numa,3*numa+1)
    hold on
    subplot(numd,numa,3*numa+2)
    hold on
    subplot(numd,numa,3*numa+3)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,53))
    title('I_{CaB}')

    subplot(numd,numa,4*numa+1)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,55))
    title('I_{TRPN}')
    subplot(numd,numa,4*numa+2)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,57))
    title('beta_{fluo}')
    subplot(numd,numa,4*numa+3)
    hold on
    plot(VOI(start:end),ALGEBRAIC(start:end,56))
    title('beta_{CMDN}')
end
end

function [algebraicVariableCount] = getAlgebraicVariableCount()
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".
    algebraicVariableCount =66;
end
% There are a total of 7 entries in each of the rate and state variable arrays.
% There are a total of 58 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over
    tspan = [0, 30*CONSTANTS(:,7)-1];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
%     [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
%     figure();
%     plot(VOI, STATES);
%     xlabel(LEGEND_VOI);
%     l = legend(LEGEND_STATES);
%     set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (ms)');
    LEGEND_CONSTANTS(:,1) = strpad('V_myo in component cell_geometry (um3)');
    LEGEND_CONSTANTS(:,2) = strpad('V_SR in component cell_geometry (um3)');
    LEGEND_CONSTANTS(:,3) = strpad('A_cap in component cell_geometry (um2)');
    LEGEND_ALGEBRAIC(:,1) = strpad('V in component membrane (mV)');
    LEGEND_CONSTANTS(:,4) = strpad('R in component membrane (mJ_per_mole_K)');
    LEGEND_CONSTANTS(:,5) = strpad('T in component membrane (kelvin)');
    LEGEND_CONSTANTS(:,6) = strpad('F in component membrane (C_per_mole)');
    LEGEND_ALGEBRAIC(:,4) = strpad('FVRT in component membrane (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('FVRT_Ca in component membrane (dimensionless)');
    LEGEND_CONSTANTS(:,7) = strpad('beat_period in component membrane (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('g_D in component CaRU (um3_per_ms)');
    LEGEND_CONSTANTS(:,9) = strpad('J_R in component CaRU (um3_per_ms)');
    LEGEND_CONSTANTS(:,10) = strpad('J_L in component CaRU (um3_per_ms)');
    LEGEND_CONSTANTS(:,11) = strpad('N in component CaRU (dimensionless)');
    LEGEND_STATES(:,1) = strpad('Ca_i in component intracell (mM)');
    LEGEND_CONSTANTS(:,12) = strpad('Ca_o in component extracellular_ion_concentrations (mM)');
    LEGEND_STATES(:,2) = strpad('Ca_SR in component intracell (mM)');
    LEGEND_ALGEBRAIC(:,35) = strpad('I_RyR in component RyR_current (mM_per_ms)');
    LEGEND_ALGEBRAIC(:,44) = strpad('I_LCC in component LCC_current (mM_per_ms)');
    LEGEND_ALGEBRAIC(:,16) = strpad('C_oc in component DS_Calcium_Concentrations (mM)');
    LEGEND_ALGEBRAIC(:,14) = strpad('C_co in component DS_Calcium_Concentrations (mM)');
    LEGEND_CONSTANTS(:,13) = strpad('V_L in component CaRU_Transitions (mV)');
    LEGEND_CONSTANTS(:,14) = strpad('del_VL in component CaRU_Transitions (mV)');
    LEGEND_CONSTANTS(:,15) = strpad('phi_L in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,16) = strpad('t_L in component CaRU_Transitions (ms)');
    LEGEND_CONSTANTS(:,17) = strpad('tau_L in component CaRU_Transitions (ms)');
    LEGEND_CONSTANTS(:,56) = strpad('t_R in component CaRU_Transitions (ms)');
    LEGEND_CONSTANTS(:,18) = strpad('tau_R in component CaRU_Transitions (ms)');
    LEGEND_CONSTANTS(:,19) = strpad('phi_R in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,20) = strpad('theta_R in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,21) = strpad('K_RyR in component CaRU_Transitions (mM)');
    LEGEND_CONSTANTS(:,22) = strpad('K_L in component CaRU_Transitions (mM)');
    LEGEND_CONSTANTS(:,23) = strpad('a in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,24) = strpad('b in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,25) = strpad('c in component CaRU_Transitions (dimensionless)');
    LEGEND_CONSTANTS(:,26) = strpad('d in component CaRU_Transitions (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('expVL in component CaRU_Transitions (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('alpha_p in component CaRU_Transitions (per_ms)');
    LEGEND_CONSTANTS(:,57) = strpad('alpha_m in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,17) = strpad('beta_poc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,9) = strpad('beta_pcc in component CaRU_Transitions (per_ms)');
    LEGEND_CONSTANTS(:,58) = strpad('beta_m in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,15) = strpad('epsilon_pco in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,10) = strpad('epsilon_pcc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,11) = strpad('epsilon_m in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,18) = strpad('mu_poc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,12) = strpad('mu_pcc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,19) = strpad('mu_moc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,13) = strpad('mu_mcc in component CaRU_Transitions (per_ms)');
    LEGEND_ALGEBRAIC(:,2) = strpad('C_cc in component DS_Calcium_Concentrations (mM)');
    LEGEND_ALGEBRAIC(:,7) = strpad('C_oo in component DS_Calcium_Concentrations (mM)');
    LEGEND_ALGEBRAIC(:,23) = strpad('J_Loo in component LCC_and_RyR_fluxes (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,22) = strpad('J_Loc in component LCC_and_RyR_fluxes (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,20) = strpad('J_Rco in component LCC_and_RyR_fluxes (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,21) = strpad('J_Roo in component LCC_and_RyR_fluxes (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,24) = strpad('denom in component CaRU_states (per_ms3)');
    LEGEND_ALGEBRAIC(:,25) = strpad('y_oc in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,26) = strpad('y_co in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,27) = strpad('y_oo in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('y_cc in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,30) = strpad('y_ci in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,33) = strpad('y_oi in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,36) = strpad('y_ic in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,39) = strpad('y_io in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,42) = strpad('y_ii in component CaRU_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,31) = strpad('r_1 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,34) = strpad('r_2 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,37) = strpad('r_3 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,40) = strpad('r_4 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,43) = strpad('r_5 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,45) = strpad('r_6 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,47) = strpad('r_7 in component CaRU_reduced_states (per_ms)');
    LEGEND_ALGEBRAIC(:,49) = strpad('r_8 in component CaRU_reduced_states (per_ms)');
    LEGEND_STATES(:,3) = strpad('z_1 in component CaRU_reduced_states (dimensionless)');
    LEGEND_STATES(:,4) = strpad('z_2 in component CaRU_reduced_states (dimensionless)');
    LEGEND_STATES(:,5) = strpad('z_3 in component CaRU_reduced_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,51) = strpad('z_4 in component CaRU_reduced_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,29) = strpad('J_R1 in component RyR_current (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,32) = strpad('J_R3 in component RyR_current (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,38) = strpad('J_L1 in component LCC_current (um3_mM_per_ms)');
    LEGEND_ALGEBRAIC(:,41) = strpad('J_L2 in component LCC_current (um3_mM_per_ms)');
    LEGEND_CONSTANTS(:,27) = strpad('K_mNa in component Na_Ca_Exchanger (mM)');
    LEGEND_CONSTANTS(:,28) = strpad('K_mCa in component Na_Ca_Exchanger (mM)');
    LEGEND_CONSTANTS(:,29) = strpad('eta in component Na_Ca_Exchanger (dimensionless)');
    LEGEND_CONSTANTS(:,30) = strpad('k_sat in component Na_Ca_Exchanger (dimensionless)');
    LEGEND_CONSTANTS(:,31) = strpad('g_NCX in component Na_Ca_Exchanger (mM_per_ms)');
    LEGEND_CONSTANTS(:,32) = strpad('Na_i in component intracell (mM)');
    LEGEND_CONSTANTS(:,33) = strpad('Na_o in component extracellular_ion_concentrations (mM)');
    LEGEND_ALGEBRAIC(:,46) = strpad('I_NaCa in component Na_Ca_Exchanger (mM_per_ms)');
    LEGEND_CONSTANTS(:,34) = strpad('g_SERCA in component SERCA (mM_per_ms)');
    LEGEND_CONSTANTS(:,35) = strpad('K_SERCA in component SERCA (mM)');
    LEGEND_ALGEBRAIC(:,48) = strpad('I_SERCA in component SERCA (mM_per_ms)');
    LEGEND_CONSTANTS(:,36) = strpad('g_pCa in component Sarcolemmal_Ca_pump (mM_per_ms)');
    LEGEND_CONSTANTS(:,37) = strpad('K_pCa in component Sarcolemmal_Ca_pump (mM)');
    LEGEND_ALGEBRAIC(:,50) = strpad('I_pCa in component Sarcolemmal_Ca_pump (mM_per_ms)');
    LEGEND_ALGEBRAIC(:,52) = strpad('E_Ca in component Background_Ca_current (mV)');
    LEGEND_CONSTANTS(:,38) = strpad('g_CaB in component Background_Ca_current (mM_per_mV_ms)');
    LEGEND_ALGEBRAIC(:,53) = strpad('I_CaB in component Background_Ca_current (mM_per_ms)');
    LEGEND_CONSTANTS(:,39) = strpad('g_SRl in component SR_Ca_leak_current (per_ms)');
    LEGEND_ALGEBRAIC(:,54) = strpad('I_SR in component SR_Ca_leak_current (mM_per_ms)');
    LEGEND_CONSTANTS(:,40) = strpad('k_m_TRPN in component troponin_Ca_buffer (per_ms)');
    LEGEND_CONSTANTS(:,41) = strpad('k_p_TRPN in component troponin_Ca_buffer (per_mM_ms)');
    LEGEND_CONSTANTS(:,42) = strpad('B_TRPN in component troponin_Ca_buffer (mM)');
    LEGEND_STATES(:,6) = strpad('TRPN in component intracell (mM)');
    LEGEND_ALGEBRAIC(:,55) = strpad('I_TRPN in component troponin_Ca_buffer (mM_per_ms)');
    LEGEND_CONSTANTS(:,43) = strpad('k_CMDN in component calmodulin_Ca_buffer (mM)');
    LEGEND_CONSTANTS(:,44) = strpad('B_CMDN in component calmodulin_Ca_buffer (mM)');
    LEGEND_ALGEBRAIC(:,56) = strpad('beta_CMDN in component calmodulin_Ca_buffer (dimensionless)');
    LEGEND_CONSTANTS(:,45) = strpad('k_fluo in component fluo4AM_buffer (mM)');
    LEGEND_CONSTANTS(:,46) = strpad('B_fluo in component fluo4AM_buffer (mM)');
    LEGEND_ALGEBRAIC(:,57) = strpad('beta_fluo in component fluo4AM_buffer (dimensionless)');
    LEGEND_ALGEBRAIC(:,3) = strpad('CaSR_plot in component intracell (mM)');
    LEGEND_ALGEBRAIC(:,65) = strpad('open_probability in component open_probability (dimensionless)');
    LEGEND_ALGEBRAIC(:,66) = strpad('I_IP3R in component intracell (mM_per_ms)');
    LEGEND_CONSTANTS(:,47) = strpad('k_f in component intracell (um3_per_ms)');
    LEGEND_ALGEBRAIC(:,58) = strpad('I_other in component intracell (mM)');
    LEGEND_CONSTANTS(:,48) = strpad('N_IP3R in component intracell (dimensionless)');
    LEGEND_ALGEBRAIC(:,64) = strpad('beta in component alphabeta (dimensionless)');
    LEGEND_ALGEBRAIC(:,63) = strpad('alpha in component alphabeta (dimensionless)');
    LEGEND_CONSTANTS(:,49) = strpad('k_beta in component open_probability (dimensionless)');
    LEGEND_ALGEBRAIC(:,62) = strpad('B in component alphabeta (dimensionless)');
    LEGEND_ALGEBRAIC(:,59) = strpad('m in component alphabeta (dimensionless)');
    LEGEND_STATES(:,7) = strpad('h in component alphabeta (dimensionless)');
    LEGEND_ALGEBRAIC(:,60) = strpad('h_alpha in component alphabeta (dimensionless)');
    LEGEND_CONSTANTS(:,50) = strpad('K_c in component alphabeta (mM)');
    LEGEND_CONSTANTS(:,51) = strpad('K_h in component alphabeta (mM)');
    LEGEND_CONSTANTS(:,52) = strpad('K_p in component alphabeta (mM)');
    LEGEND_ALGEBRAIC(:,61) = strpad('IP_3 in component alphabeta (mM)');
    LEGEND_CONSTANTS(:,53) = strpad('t_max in component alphabeta (per_ms)');
    LEGEND_CONSTANTS(:,54) = strpad('K_t in component alphabeta (mM)');
    LEGEND_CONSTANTS(:,55) = strpad('t_p in component alphabeta (per_ms)');
    LEGEND_RATES(:,3) = strpad('d/dt z_1 in component CaRU_reduced_states (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt z_2 in component CaRU_reduced_states (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt z_3 in component CaRU_reduced_states (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt TRPN in component intracell (mM)');
    LEGEND_RATES(:,1) = strpad('d/dt Ca_i in component intracell (mM)');
    LEGEND_RATES(:,2) = strpad('d/dt Ca_SR in component intracell (mM)');
    LEGEND_RATES(:,7) = strpad('d/dt h in component alphabeta (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    global var1;global var3; global var4; global var5;global var6;global initStates;
    if isempty(initStates)
        STATES(:,:)=[0.000102316812186319,0.683347094971423,0.988494715189347,0.00875264158667971,0.00272874220458962,0.0634814538696518,0.240238948161970];
    else
        STATES(:,:)=initStates;
    end
    CONSTANTS(:,1) = 25.84e3;
    CONSTANTS(:,2) = 2.098e3;
    CONSTANTS(:,3) = 1.534e4;
    CONSTANTS(:,4) = 8314.5;
    CONSTANTS(:,5) = 295;
    CONSTANTS(:,6) = 96487;
    CONSTANTS(:,7) = 1000;
    CONSTANTS(:,8) = 0.065;
    CONSTANTS(:,9) = 0.02;
    CONSTANTS(:,10) = 9.13e-4;
    CONSTANTS(:,11) = 50000;
    CONSTANTS(:,12) = 1;
    CONSTANTS(:,13) = -2;
    CONSTANTS(:,14) = 7;
    CONSTANTS(:,15) = 2.35;
    CONSTANTS(:,16) = 1;
    CONSTANTS(:,17) = 650;
    CONSTANTS(:,18) = 2.43;
    CONSTANTS(:,19) = 0.05;
    CONSTANTS(:,20) = 0.012;
    CONSTANTS(:,21) = 41e-3;
    CONSTANTS(:,22) = 0.22e-3;
    CONSTANTS(:,23) = 0.0625;
    CONSTANTS(:,24) = 14;
    CONSTANTS(:,25) = 0.01;
    CONSTANTS(:,26) = 100;
    CONSTANTS(:,27) = 87.5;
    CONSTANTS(:,28) = 1.38;
    CONSTANTS(:,29) = 0.35;
    CONSTANTS(:,30) = 0.1;
    CONSTANTS(:,31) = 38.5e-3;
    CONSTANTS(:,32) = 10;
    CONSTANTS(:,33) = 140;
    CONSTANTS(:,34) = 0.45e-3;
    CONSTANTS(:,35) = 0.5e-3;
    CONSTANTS(:,36) = 0.0035e-3;
    CONSTANTS(:,37) = 0.5e-3;
    CONSTANTS(:,38) = 2.6875e-8;
    CONSTANTS(:,39) = 1.8951e-5;
    CONSTANTS(:,40) = 0.04;
    CONSTANTS(:,41) = 0.04e3;
    CONSTANTS(:,42) = 70e-3;
    CONSTANTS(:,43) = 2.382e-3;
    CONSTANTS(:,44) = 50e-3;
    CONSTANTS(:,45) = 1;
    CONSTANTS(:,46) = 1e-3;
    CONSTANTS(:,47) = 3e-2*var1;
    CONSTANTS(:,48) = 1000;
    CONSTANTS(:,49) = 0.324;
    CONSTANTS(:,50) = 2e-4*var4;%K_c
    CONSTANTS(:,51) = 8e-5*var5; %K_h
    CONSTANTS(:,52) = 2e-4; 
    CONSTANTS(:,53) = 1*var3;
    CONSTANTS(:,54) = 1e-4*var6;
    CONSTANTS(:,55) = 2.7e-5;
    CONSTANTS(:,56) =  1.17000.*CONSTANTS(:,16);
    CONSTANTS(:,57) = CONSTANTS(:,15)./CONSTANTS(:,16);
    CONSTANTS(:,58) = CONSTANTS(:,19)./CONSTANTS(:,56);
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount; global var1;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,1) = piecewise({ rem(VOI, CONSTANTS(:,7))>=0.00000& rem(VOI, CONSTANTS(:,7))<=80.000, 35-1.4375*rem(VOI, CONSTANTS(:,7))},  - 80.0000);
    ALGEBRAIC(:,4) = ( CONSTANTS(:,6).*ALGEBRAIC(:,1))./( CONSTANTS(:,4).*CONSTANTS(:,5));
    ALGEBRAIC(:,5) =  2.00000.*ALGEBRAIC(:,4);
    ALGEBRAIC(:,16) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-09, (STATES(:,1)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12).*ALGEBRAIC(:,5).*exp( - ALGEBRAIC(:,5)))./(1.00000 - exp( - ALGEBRAIC(:,5))))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, (STATES(:,1)+ (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12))./(1.00000+CONSTANTS(:,10)./CONSTANTS(:,8)));
    ALGEBRAIC(:,18) = (power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000))./( CONSTANTS(:,18).*(power(ALGEBRAIC(:,16), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,12) = (power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000))./( CONSTANTS(:,18).*(power(STATES(:,1), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,6) = exp((ALGEBRAIC(:,1) - CONSTANTS(:,13))./CONSTANTS(:,14));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,6)./( CONSTANTS(:,16).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,9) = power(STATES(:,1), 2.00000)./( CONSTANTS(:,56).*(power(STATES(:,1), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,17) = power(ALGEBRAIC(:,16), 2.00000)./( CONSTANTS(:,56).*(power(ALGEBRAIC(:,16), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,24) =  (ALGEBRAIC(:,8)+CONSTANTS(:,57)).*( (CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,17)).*(CONSTANTS(:,58)+ALGEBRAIC(:,9))+ ALGEBRAIC(:,8).*(CONSTANTS(:,58)+ALGEBRAIC(:,17)));
    ALGEBRAIC(:,25) = ( ALGEBRAIC(:,8).*CONSTANTS(:,58).*(ALGEBRAIC(:,8)+CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,9)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,28) = ( CONSTANTS(:,57).*CONSTANTS(:,58).*(CONSTANTS(:,57)+ALGEBRAIC(:,8)+CONSTANTS(:,58)+ALGEBRAIC(:,17)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,31) =  ALGEBRAIC(:,25).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,28).*ALGEBRAIC(:,12);
    ALGEBRAIC(:,19) = ( CONSTANTS(:,20).*CONSTANTS(:,26).*(power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)))./( CONSTANTS(:,18).*( CONSTANTS(:,26).*power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,13) = ( CONSTANTS(:,20).*CONSTANTS(:,26).*(power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)))./( CONSTANTS(:,18).*( CONSTANTS(:,26).*power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,34) = ( ALGEBRAIC(:,8).*ALGEBRAIC(:,19)+ CONSTANTS(:,57).*ALGEBRAIC(:,13))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,14) = (STATES(:,1)+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*STATES(:,2))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8));
    ALGEBRAIC(:,15) = ( ALGEBRAIC(:,14).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*CONSTANTS(:,22).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,10) = ( STATES(:,1).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*CONSTANTS(:,22).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,26) = ( CONSTANTS(:,57).*( ALGEBRAIC(:,9).*(CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,17))+ ALGEBRAIC(:,17).*ALGEBRAIC(:,8)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,43) =  ALGEBRAIC(:,26).*ALGEBRAIC(:,15)+ ALGEBRAIC(:,28).*ALGEBRAIC(:,10);
    ALGEBRAIC(:,11) = ( CONSTANTS(:,24).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*( CONSTANTS(:,24).*ALGEBRAIC(:,6)+CONSTANTS(:,23)));
    ALGEBRAIC(:,45) = ALGEBRAIC(:,11);
    RATES(:,3) =   - (ALGEBRAIC(:,31)+ALGEBRAIC(:,43)).*STATES(:,3)+ ALGEBRAIC(:,34).*STATES(:,4)+ ALGEBRAIC(:,45).*STATES(:,5);
    ALGEBRAIC(:,47) = ( CONSTANTS(:,57).*ALGEBRAIC(:,10))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,49) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,51) = ((1.00000 - STATES(:,3)) - STATES(:,4)) - STATES(:,5);
    RATES(:,4) = ( ALGEBRAIC(:,31).*STATES(:,3) -  (ALGEBRAIC(:,34)+ALGEBRAIC(:,47)).*STATES(:,4))+ ALGEBRAIC(:,49).*ALGEBRAIC(:,51);
    ALGEBRAIC(:,37) = ( CONSTANTS(:,58).*ALGEBRAIC(:,12))./(CONSTANTS(:,58)+ALGEBRAIC(:,9));
    ALGEBRAIC(:,40) = ALGEBRAIC(:,13);
    RATES(:,5) = ( ALGEBRAIC(:,43).*STATES(:,3) -  (ALGEBRAIC(:,45)+ALGEBRAIC(:,37)).*STATES(:,5))+ ALGEBRAIC(:,40).*ALGEBRAIC(:,51);
    ALGEBRAIC(:,55) =  CONSTANTS(:,40).*(CONSTANTS(:,42) - STATES(:,6)) -  CONSTANTS(:,41).*STATES(:,6).*STATES(:,1);
    RATES(:,6) = ALGEBRAIC(:,55);
    ALGEBRAIC(:,60) = power(CONSTANTS(:,51), 4.00000)./(power(CONSTANTS(:,51), 4.00000)+power(STATES(:,1), 4.00000));
    RATES(:,7) = (var1~=0)*( (ALGEBRAIC(:,60) - STATES(:,7)).*(power(CONSTANTS(:,54), 4.00000)+power(STATES(:,1), 4.00000)))./( CONSTANTS(:,53).*power(CONSTANTS(:,54), 4.00000));
    ALGEBRAIC(:,20) = ( CONSTANTS(:,9).*(STATES(:,2) - STATES(:,1)))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8));
    ALGEBRAIC(:,21) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( CONSTANTS(:,9).*((STATES(:,2) - STATES(:,1))+ (( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*(STATES(:,2) -  CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, ( CONSTANTS(:,9).*((STATES(:,2) - STATES(:,1))+ (( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(STATES(:,2) -  CONSTANTS(:,12).*exp( - 1.00000e-05))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,27) = ( ALGEBRAIC(:,8).*( ALGEBRAIC(:,17).*(ALGEBRAIC(:,8)+CONSTANTS(:,58)+ALGEBRAIC(:,9))+ ALGEBRAIC(:,9).*CONSTANTS(:,57)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,29) =  ALGEBRAIC(:,27).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,20).*ALGEBRAIC(:,26);
    ALGEBRAIC(:,32) = ( ALGEBRAIC(:,20).*ALGEBRAIC(:,9))./(CONSTANTS(:,58)+ALGEBRAIC(:,9));
    ALGEBRAIC(:,35) = ( ( STATES(:,3).*ALGEBRAIC(:,29)+ STATES(:,5).*ALGEBRAIC(:,32)).*CONSTANTS(:,11))./CONSTANTS(:,1);
    ALGEBRAIC(:,48) = ( CONSTANTS(:,34).*power(STATES(:,1), 2.00000))./(power(CONSTANTS(:,35), 2.00000)+power(STATES(:,1), 2.00000));
    ALGEBRAIC(:,56) = power(1.00000+( CONSTANTS(:,43).*CONSTANTS(:,44))./power(CONSTANTS(:,43)+STATES(:,1), 2.00000),  - 1.00000);
    ALGEBRAIC(:,57) = power(1.00000+( CONSTANTS(:,45).*CONSTANTS(:,46))./power(CONSTANTS(:,45)+STATES(:,1), 2.00000),  - 1.00000);
    ALGEBRAIC(:,61) = 1.00000e-02;%*var1;%piecewise({VOI<3300.00, 0.00000 }, 1.00000e-02);
    ALGEBRAIC(:,62) = power(ALGEBRAIC(:,61), 2.00000)./(power(CONSTANTS(:,52), 2.00000)+power(ALGEBRAIC(:,61), 2.00000));
    ALGEBRAIC(:,59) = power(STATES(:,1), 4.00000)./(power(CONSTANTS(:,50), 4.00000)+power(STATES(:,1), 4.00000));
    ALGEBRAIC(:,64) =  ALGEBRAIC(:,62).*ALGEBRAIC(:,59).*STATES(:,7);
    ALGEBRAIC(:,63) =  (1.00000 - ALGEBRAIC(:,62)).*(1.00000 -  ALGEBRAIC(:,59).*ALGEBRAIC(:,60));
    ALGEBRAIC(:,65) = ALGEBRAIC(:,64)./(ALGEBRAIC(:,64)+ CONSTANTS(:,49).*(ALGEBRAIC(:,64)+ALGEBRAIC(:,63)));
    ALGEBRAIC(:,66) = ( CONSTANTS(:,47).*CONSTANTS(:,48).*ALGEBRAIC(:,65).*(STATES(:,2) - STATES(:,1)))./CONSTANTS(:,1);
    ALGEBRAIC(:,23) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( (( CONSTANTS(:,10).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*(( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,1))+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,2))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp(ALGEBRAIC(:,5)))) }, ( (( CONSTANTS(:,10).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,1))+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,2))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,22) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( (( CONSTANTS(:,10).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,1)))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, ( (( CONSTANTS(:,10).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,1)))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,38) =  ALGEBRAIC(:,23).*ALGEBRAIC(:,27)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,25);
    ALGEBRAIC(:,41) = ( ALGEBRAIC(:,22).*ALGEBRAIC(:,8))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,44) = ( ( STATES(:,3).*ALGEBRAIC(:,38)+ STATES(:,4).*ALGEBRAIC(:,41)).*CONSTANTS(:,11))./CONSTANTS(:,1);
    ALGEBRAIC(:,46) = ( CONSTANTS(:,31).*( exp( CONSTANTS(:,29).*ALGEBRAIC(:,4)).*power(CONSTANTS(:,32), 3.00000).*CONSTANTS(:,12) -  exp( (CONSTANTS(:,29) - 1.00000).*ALGEBRAIC(:,4)).*power(CONSTANTS(:,33), 3.00000).*STATES(:,1)))./( (power(CONSTANTS(:,33), 3.00000)+power(CONSTANTS(:,27), 3.00000)).*(CONSTANTS(:,12)+CONSTANTS(:,28)).*(1.00000+ CONSTANTS(:,30).*exp( (CONSTANTS(:,29) - 1.00000).*ALGEBRAIC(:,4))));
    ALGEBRAIC(:,50) = ( CONSTANTS(:,36).*STATES(:,1))./(CONSTANTS(:,37)+STATES(:,1));
    ALGEBRAIC(:,52) =  (( CONSTANTS(:,4).*CONSTANTS(:,5))./( 2.00000.*CONSTANTS(:,6))).*log(CONSTANTS(:,12)./STATES(:,1));
    ALGEBRAIC(:,53) =  CONSTANTS(:,38).*(ALGEBRAIC(:,52) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,54) =  CONSTANTS(:,39).*(STATES(:,2) - STATES(:,1));
    ALGEBRAIC(:,58) = ((ALGEBRAIC(:,44)+ALGEBRAIC(:,54)+ALGEBRAIC(:,46)) - ALGEBRAIC(:,50))+ALGEBRAIC(:,53)+ALGEBRAIC(:,55);
    RATES(:,1) =  ALGEBRAIC(:,57).*ALGEBRAIC(:,56).*((ALGEBRAIC(:,35) - ALGEBRAIC(:,48))+ALGEBRAIC(:,66)+ALGEBRAIC(:,58));
    RATES(:,2) =  (CONSTANTS(:,1)./CONSTANTS(:,2)).*((( - ALGEBRAIC(:,35)+ALGEBRAIC(:,48)) - ALGEBRAIC(:,54)) - ALGEBRAIC(:,66));
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    global var1;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,1) = piecewise({ rem(VOI, CONSTANTS(:,7))>=0.00000& rem(VOI, CONSTANTS(:,7))<=80.000, 35-1.4375*rem(VOI, CONSTANTS(:,7))},  - 80.0000);
    ALGEBRAIC(:,4) = ( CONSTANTS(:,6).*ALGEBRAIC(:,1))./( CONSTANTS(:,4).*CONSTANTS(:,5));
    ALGEBRAIC(:,5) =  2.00000.*ALGEBRAIC(:,4);
    ALGEBRAIC(:,16) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-09, (STATES(:,1)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12).*ALGEBRAIC(:,5).*exp( - ALGEBRAIC(:,5)))./(1.00000 - exp( - ALGEBRAIC(:,5))))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, (STATES(:,1)+ (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12))./(1.00000+CONSTANTS(:,10)./CONSTANTS(:,8)));
    ALGEBRAIC(:,18) = (power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000))./( CONSTANTS(:,18).*(power(ALGEBRAIC(:,16), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,12) = (power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000))./( CONSTANTS(:,18).*(power(STATES(:,1), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,6) = exp((ALGEBRAIC(:,1) - CONSTANTS(:,13))./CONSTANTS(:,14));
    ALGEBRAIC(:,8) = ALGEBRAIC(:,6)./( CONSTANTS(:,16).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,9) = power(STATES(:,1), 2.00000)./( CONSTANTS(:,56).*(power(STATES(:,1), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,17) = power(ALGEBRAIC(:,16), 2.00000)./( CONSTANTS(:,56).*(power(ALGEBRAIC(:,16), 2.00000)+power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,24) =  (ALGEBRAIC(:,8)+CONSTANTS(:,57)).*( (CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,17)).*(CONSTANTS(:,58)+ALGEBRAIC(:,9))+ ALGEBRAIC(:,8).*(CONSTANTS(:,58)+ALGEBRAIC(:,17)));
    ALGEBRAIC(:,25) = ( ALGEBRAIC(:,8).*CONSTANTS(:,58).*(ALGEBRAIC(:,8)+CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,9)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,28) = ( CONSTANTS(:,57).*CONSTANTS(:,58).*(CONSTANTS(:,57)+ALGEBRAIC(:,8)+CONSTANTS(:,58)+ALGEBRAIC(:,17)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,31) =  ALGEBRAIC(:,25).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,28).*ALGEBRAIC(:,12);
    ALGEBRAIC(:,19) = ( CONSTANTS(:,20).*CONSTANTS(:,26).*(power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)))./( CONSTANTS(:,18).*( CONSTANTS(:,26).*power(ALGEBRAIC(:,16), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,13) = ( CONSTANTS(:,20).*CONSTANTS(:,26).*(power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)))./( CONSTANTS(:,18).*( CONSTANTS(:,26).*power(STATES(:,1), 2.00000)+ CONSTANTS(:,25).*power(CONSTANTS(:,21), 2.00000)));
    ALGEBRAIC(:,34) = ( ALGEBRAIC(:,8).*ALGEBRAIC(:,19)+ CONSTANTS(:,57).*ALGEBRAIC(:,13))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,14) = (STATES(:,1)+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*STATES(:,2))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8));
    ALGEBRAIC(:,15) = ( ALGEBRAIC(:,14).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*CONSTANTS(:,22).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,10) = ( STATES(:,1).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*CONSTANTS(:,22).*(ALGEBRAIC(:,6)+1.00000));
    ALGEBRAIC(:,26) = ( CONSTANTS(:,57).*( ALGEBRAIC(:,9).*(CONSTANTS(:,57)+CONSTANTS(:,58)+ALGEBRAIC(:,17))+ ALGEBRAIC(:,17).*ALGEBRAIC(:,8)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,43) =  ALGEBRAIC(:,26).*ALGEBRAIC(:,15)+ ALGEBRAIC(:,28).*ALGEBRAIC(:,10);
    ALGEBRAIC(:,11) = ( CONSTANTS(:,24).*(ALGEBRAIC(:,6)+CONSTANTS(:,23)))./( CONSTANTS(:,17).*( CONSTANTS(:,24).*ALGEBRAIC(:,6)+CONSTANTS(:,23)));
    ALGEBRAIC(:,45) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,47) = ( CONSTANTS(:,57).*ALGEBRAIC(:,10))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,49) = ALGEBRAIC(:,11);
    ALGEBRAIC(:,51) = ((1.00000 - STATES(:,3)) - STATES(:,4)) - STATES(:,5);
    ALGEBRAIC(:,37) = ( CONSTANTS(:,58).*ALGEBRAIC(:,12))./(CONSTANTS(:,58)+ALGEBRAIC(:,9));
    ALGEBRAIC(:,40) = ALGEBRAIC(:,13);
    ALGEBRAIC(:,55) =  CONSTANTS(:,40).*(CONSTANTS(:,42) - STATES(:,6)) -  CONSTANTS(:,41).*STATES(:,6).*STATES(:,1);
    ALGEBRAIC(:,60) = power(CONSTANTS(:,51), 4.00000)./(power(CONSTANTS(:,51), 4.00000)+power(STATES(:,1), 4.00000));
    ALGEBRAIC(:,20) = ( CONSTANTS(:,9).*(STATES(:,2) - STATES(:,1)))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8));
    ALGEBRAIC(:,21) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( CONSTANTS(:,9).*((STATES(:,2) - STATES(:,1))+ (( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*(STATES(:,2) -  CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, ( CONSTANTS(:,9).*((STATES(:,2) - STATES(:,1))+ (( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(STATES(:,2) -  CONSTANTS(:,12).*exp( - 1.00000e-05))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,27) = ( ALGEBRAIC(:,8).*( ALGEBRAIC(:,17).*(ALGEBRAIC(:,8)+CONSTANTS(:,58)+ALGEBRAIC(:,9))+ ALGEBRAIC(:,9).*CONSTANTS(:,57)))./ALGEBRAIC(:,24);
    ALGEBRAIC(:,29) =  ALGEBRAIC(:,27).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,20).*ALGEBRAIC(:,26);
    ALGEBRAIC(:,32) = ( ALGEBRAIC(:,20).*ALGEBRAIC(:,9))./(CONSTANTS(:,58)+ALGEBRAIC(:,9));
    ALGEBRAIC(:,35) = ( ( STATES(:,3).*ALGEBRAIC(:,29)+ STATES(:,5).*ALGEBRAIC(:,32)).*CONSTANTS(:,11))./CONSTANTS(:,1);
    ALGEBRAIC(:,48) = ( CONSTANTS(:,34).*power(STATES(:,1), 2.00000))./(power(CONSTANTS(:,35), 2.00000)+power(STATES(:,1), 2.00000));
    ALGEBRAIC(:,56) = power(1.00000+( CONSTANTS(:,43).*CONSTANTS(:,44))./power(CONSTANTS(:,43)+STATES(:,1), 2.00000),  - 1.00000);
    ALGEBRAIC(:,57) = power(1.00000+( CONSTANTS(:,45).*CONSTANTS(:,46))./power(CONSTANTS(:,45)+STATES(:,1), 2.00000),  - 1.00000);
    ALGEBRAIC(:,61) = 1.00000e-02;%*var1;%piecewise({VOI<3300.00, 0.00000 }, 1.00000e-02);
    ALGEBRAIC(:,62) = power(ALGEBRAIC(:,61), 2.00000)./(power(CONSTANTS(:,52), 2.00000)+power(ALGEBRAIC(:,61), 2.00000));
    ALGEBRAIC(:,59) = power(STATES(:,1), 4.00000)./(power(CONSTANTS(:,50), 4.00000)+power(STATES(:,1), 4.00000));
    ALGEBRAIC(:,64) =  ALGEBRAIC(:,62).*ALGEBRAIC(:,59).*STATES(:,7);
    ALGEBRAIC(:,63) =  (1.00000 - ALGEBRAIC(:,62)).*(1.00000 -  ALGEBRAIC(:,59).*ALGEBRAIC(:,60));
    ALGEBRAIC(:,65) = ALGEBRAIC(:,64)./(ALGEBRAIC(:,64)+ CONSTANTS(:,49).*(ALGEBRAIC(:,64)+ALGEBRAIC(:,63)));
    ALGEBRAIC(:,66) = ( CONSTANTS(:,47).*CONSTANTS(:,48).*ALGEBRAIC(:,65).*(STATES(:,2) - STATES(:,1)))./CONSTANTS(:,1);
    ALGEBRAIC(:,23) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( (( CONSTANTS(:,10).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*(( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,1))+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,2))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp(ALGEBRAIC(:,5)))) }, ( (( CONSTANTS(:,10).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*(( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,1))+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,2))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,22) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-05, ( (( CONSTANTS(:,10).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))).*( CONSTANTS(:,12).*exp( - ALGEBRAIC(:,5)) - STATES(:,1)))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, ( (( CONSTANTS(:,10).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))).*( CONSTANTS(:,12).*exp( - 1.00000e-05) - STATES(:,1)))./(1.00000+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*1.00000e-05)./(1.00000 - exp( - 1.00000e-05))));
    ALGEBRAIC(:,38) =  ALGEBRAIC(:,23).*ALGEBRAIC(:,27)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,25);
    ALGEBRAIC(:,41) = ( ALGEBRAIC(:,22).*ALGEBRAIC(:,8))./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,44) = ( ( STATES(:,3).*ALGEBRAIC(:,38)+ STATES(:,4).*ALGEBRAIC(:,41)).*CONSTANTS(:,11))./CONSTANTS(:,1);
    ALGEBRAIC(:,46) = ( CONSTANTS(:,31).*( exp( CONSTANTS(:,29).*ALGEBRAIC(:,4)).*power(CONSTANTS(:,32), 3.00000).*CONSTANTS(:,12) -  exp( (CONSTANTS(:,29) - 1.00000).*ALGEBRAIC(:,4)).*power(CONSTANTS(:,33), 3.00000).*STATES(:,1)))./( (power(CONSTANTS(:,33), 3.00000)+power(CONSTANTS(:,27), 3.00000)).*(CONSTANTS(:,12)+CONSTANTS(:,28)).*(1.00000+ CONSTANTS(:,30).*exp( (CONSTANTS(:,29) - 1.00000).*ALGEBRAIC(:,4))));
    ALGEBRAIC(:,50) = ( CONSTANTS(:,36).*STATES(:,1))./(CONSTANTS(:,37)+STATES(:,1));
    ALGEBRAIC(:,52) =  (( CONSTANTS(:,4).*CONSTANTS(:,5))./( 2.00000.*CONSTANTS(:,6))).*log(CONSTANTS(:,12)./STATES(:,1));
    ALGEBRAIC(:,53) =  CONSTANTS(:,38).*(ALGEBRAIC(:,52) - ALGEBRAIC(:,1));
    ALGEBRAIC(:,54) =  CONSTANTS(:,39).*(STATES(:,2) - STATES(:,1));
    ALGEBRAIC(:,58) = ((ALGEBRAIC(:,44)+ALGEBRAIC(:,54)+ALGEBRAIC(:,46)) - ALGEBRAIC(:,50))+ALGEBRAIC(:,53)+ALGEBRAIC(:,55);
    ALGEBRAIC(:,2) = STATES(:,1);
    ALGEBRAIC(:,3) = ( STATES(:,2).*CONSTANTS(:,2))./CONSTANTS(:,1);
    ALGEBRAIC(:,7) = piecewise({abs(ALGEBRAIC(:,5))>1.00000e-09, (STATES(:,1)+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*STATES(:,2)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12).*ALGEBRAIC(:,5).*exp( - ALGEBRAIC(:,5)))./(1.00000 - exp( - ALGEBRAIC(:,5))))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+( (CONSTANTS(:,10)./CONSTANTS(:,8)).*ALGEBRAIC(:,5))./(1.00000 - exp( - ALGEBRAIC(:,5)))) }, (STATES(:,1)+ (CONSTANTS(:,9)./CONSTANTS(:,8)).*STATES(:,2)+ (CONSTANTS(:,10)./CONSTANTS(:,8)).*CONSTANTS(:,12))./(1.00000+CONSTANTS(:,9)./CONSTANTS(:,8)+CONSTANTS(:,10)./CONSTANTS(:,8)));
    ALGEBRAIC(:,30) = CONSTANTS(:,57)./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,33) = ALGEBRAIC(:,8)./(ALGEBRAIC(:,8)+CONSTANTS(:,57));
    ALGEBRAIC(:,36) = CONSTANTS(:,58)./(ALGEBRAIC(:,9)+CONSTANTS(:,58));
    ALGEBRAIC(:,39) = ALGEBRAIC(:,9)./(ALGEBRAIC(:,9)+CONSTANTS(:,58));
    ALGEBRAIC(:,42) = (((((((1.00000 - ALGEBRAIC(:,25)) - ALGEBRAIC(:,26)) - ALGEBRAIC(:,27)) - ALGEBRAIC(:,28)) - ALGEBRAIC(:,30)) - ALGEBRAIC(:,36)) - ALGEBRAIC(:,33)) - ALGEBRAIC(:,39);
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

