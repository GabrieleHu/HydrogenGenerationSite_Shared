clear all
close all
clc

%% Set up Gurobi: Adding path to optimizer

addpath('C:\gurobi1101\win64\examples\matlab')

%% Input data


% Define the main path manually, e.g. 'C:\Users\My_models'
path_main = 'C:\Users\gahu\Dropbox\PC_EMPA\Desktop\EMPA_gh\8_Coding\Green_Hydrogen\My_models\P2Xpaper_Shared\Main';
% path for input files
path_input = fullfile(path_main, 'Input');
% Define the path for functions
path_functions = fullfile(path_main, 'Functions');

% Reading an Excel file from the Input directory
excelFilePath = fullfile(path_input, 'Inputs_BoundaryLoads_2022.xlsx');
Input = readtable(excelFilePath);


%% to be modified

Input = readtable('C:\Users\gahu\Dropbox\PC_EMPA\Desktop\EMPA_gh\8_Coding\Green_Hydrogen\My_models\P2Xpaper_Shared\Main\Input\Inputs_BoundaryLoads_2022.xlsx');

%% Pre-processing

irradiance     = Input.G;                     % Hourly solar irradiance    [W/m^2]
T_amb          = Input.T_amb;                 % Hourly temperature         [°C]
cost_el        = Input.price_Eur_MWh./10^3;   % Hourly cost of electricty  [EUR/kWh]
cost_export_el = Input.Price_DayAhed./10^3;   % Hourly cost exported elec  [EUR/kWh]
heat_dem_LT    = Input.Heating; % Hourly LT heat demand                    [kW]
heat_dem_HT    = Input.DHW; % Hourly HT heat demand                        [kW]
nHours         = numel(Input.time);           % Number of hours simulated    
idxHr2ToEnd    = (2:nHours)';                 % Hours until the end
Time           = (1:nHours)';                 % Time vector
days           = nHours/24;                   % Number of days simulated
weeks          = days/7;                      % Number of weeks simulated
clear t;
linew          = 1;             
font           = 18;

% Get the efficiency and cost curve coefficients

N_bp = 2;       % Number of breakpoints on electrolyser efficiency curve: select N_bp=1 to use constant efficiency equal to nominal efficiency or N_bp=2 for linear efficiency
run efficiencies_Git.m
e_surplus = 1;  % Factor used in sensitivity analysis to scale electrolyser cost curve
run cost_curves_Git.m


%%      INPUT PARAMETERS

% general inputs

    bigM            = 1E8;     % Large number for big-M constraints
    MaxSimTime      = 120;    % Maximum time for MILP solver                  [s]
    mass_H2_day_obj = 100;     % Daily produced hydrogen                       [kg]
    H2_demand       = mass_H2_day_obj.*ones(days,1);
    cost_startup_e  = 17;     % startup penalization cost                      [EUR/h]                       
    deltat          = 3600;    % time step                                     [s]        

% Fixed technical parameters  

    HHV         = 39.39 * 3.6 * 10^3; % Hydrogen higher heating value            [kJ/kg]  
    k           = 1.4;                          % Ratio cp/cv                  [-]
    R_H2        = 4.12;                         % Universal gas constant       [kJ/kg/K] 
    T_in_H2     = 70 + 273.15;                  % H2 temperature (=T_cat=T_an) [K]  
    p_out       = 820;                          % Compressor outlet pressure   [bar]
    p_in        = 50;                           % Compressor inlet pressure    [bar]
    L_is_c      = k / (k-1) * R_H2 * T_in_H2 * ...
        (( p_out / p_in )^( ( k-1 ) /k ) -1 );  % Specific work compressor     [kJ/kg] 
    
    eff_PV      = max(eff_cell);                % (Constant) efficiency PV panels for calc
    eff_b       = 0.92;                         % Efficiency battery    
    eff_th      = 0.95;                         % Efficiency heat recovery => figure 5.3 of Tiktak master thesis
    eff_e_nom   = mm_elec(end) + qq_elec(end);  % Electrolyser efficiency at nominal power (approx)
    eff_e_c     = 0.96;                         % Electric generator efficiency compr.
    eff_m_c     = 0.98;                         % Mechanical efficiency compr.
    eff_is_c    = 0.8;                          % Isentropic efficiency compr.
    eff_c       = eff_is_c * eff_m_c * eff_e_c; % Overall efficiency compressor 0.7526

% limits equipment sizing

    S_e_min     = 100;                % Minimum size elec (feasible problem)   [kW]
    S_e_max     = 850;                % Maximum size elec                      [kW]
    C_b_max     = 600 * 3.6 * 10^3;   % Battery capacity                       [kJ]
    P_b_max     = C_b_max / 3600;     % Max C-rate battery is 1C
    mdot_H2_max = eff_e_nom * S_e_max / HHV;    % Maximum mass flow hydrogen   [kg/s]
    S_c_max     = mdot_H2_max * L_is_c / eff_c; % Maximum size compressor      [kW]
    
    P_peak_ref  = 1500;                % Reference value peak power PV        [kW]
    Area_PV_ref = P_peak_ref/eff_PV;   % Reference area PV                    [m2]
    Area_PV_min = 0;                   % Minimum area PV                      [m2]
    Area_PV_max = Area_PV_ref*2;
    

% unit prices and lifetime components

    d           = 0.04;     % Discount rate, as encouraged by EU
    ann         = d / (1 - (1 + d)^(-20));  % annuity factor calculated with plant lifetime        
    UP_PV       = 800;      % Unit price PV                                    [EUR/kW_p]
    life_PV     = 30;       % Lifetime PV                                      [years]
    maint_PV    = 0.0158;   % Annual cost maintenance PV, frac total cost
    ann_PV      = d / (1 - (1 + d)^(-life_PV));
    UP_b        = 600;      % Unit price battery                               [EUR/kWh]
    life_b      = 10;       % Lifetime battery                                 [years]
    maint_b     = 0.02;     % Annual cost maintenance battery, frac total ann cost
    ann_b       = d / (1 - (1 + d)^(-life_b));
    UP_e        = 1100;     % Unit price electrolyser                          [EUR/kW]
    cost_e_max  = aa_e(2)*S_e_max + bb_e(2);
    life_e      = 7;        % Lifetime electrolyser                            [years]
    maint_e     = 0.02;     % Annual cost maintenance electrolyser, frac total cost
    ann_e       = d / (1 - (1 + d)^(-life_e));
    UP_HP       = 576;      % Unit price heat pump                             [EUR/kW]             
    life_HP     = 20;       % Lifetime heat pump                               [years]           
    maint_HP    = 0.015;    % Annual cost maintenance HP, frac total ann cost
    ann_HP      = d / (1 - (1 + d)^(-life_HP));
    UP_HEX      = 77.79;    % Unit price heat exchanger                        [EUR/m2]         
    Fixed_HEX   = 0;        % Fixed cost for heat exchanger                    [EUR]    
    life_HEX    = 20;       % Lifetime heat exchanger                          [years]           
    maint_HEX   = 0.01;     % Annual cost maintenance HEX, frac total ann cost                   
    ann_HEX     = d / (1 - (1 + d)^(-life_HEX));
    cost_c_max  = aa_c(2)*S_c_max + bb_c(2);
    life_c      = 10;       % Lifetime compressor                              [years]           
    maint_c     = 0.08;     % Annual cost maintenance compr, frac investment cost
    ann_c       = d / (1 - (1 + d)^(-life_c));
    UP_storage  = 1644;     % Unit price hydrogen storage                      [EUR/kgH2]         
    life_storage= 10;       % Lifetime hydrogen storage                        [years]            
    % maint_storage = 0;       % NO annual cost maintenance storage
    ann_storage = d / (1 - (1 + d)^(-life_storage));
    P_refr      = 16.3;     % Refrigerator power necessary                     [kW]              
    UP_refr     = 5374;     % Unit price hydrogen refrigerator                 [EUR/kW]          
    life_refr   = 15;       % Lifetime hydrogen refrigerator                   [years]           
    maint_refr  = 0.03;     % Annual cost maintenance refrigerator, frac investment cost 
    ann_refr    = d / (1 - (1 + d)^(-life_refr));
    UP_disp     = 65000;    % Investment cost dispenser                        [EUR]             
    life_disp   = 10;       % Lifetime hydrogen dispenser                      [years]           
    maint_disp  = 0.03;     % Annual cost maintenance dispenser, frac investment cost
    ann_disp    = d / (1 - (1 + d)^(-life_disp));


% export price of heat (weighting factors from Baldini et al.)

    weight_winter = [0.432; 0.240];
    weight_mid    = [0.306; 0.107];
    weight_summer = [0.137; -0.025];
    
    cost_export_heatHT = cost_el;
    cost_export_heatLT = cost_el;
    
    winter_1 = 31*24 + 28*24 + 21*24;
    mid_1    = winter_1 + 10*24 + 30*24 + 31*24 + 21*24;
    summer   = mid_1 + 9*24 + 31*24 + 31*24 + 21*24;
    mid_2    = summer + 9*24 + 31*24 + 30*24 + 21*24;
    winter_2 = mid_2 + 10*24;
    
    cost_export_heatHT(1:winter_1)        = weight_winter(1).*cost_export_heatHT(1:winter_1);
    cost_export_heatHT(winter_1+1:mid_1)  = weight_mid(1).*cost_export_heatHT(winter_1+1:mid_1);
    cost_export_heatHT(mid_1+1:summer)    = weight_summer(1).*cost_export_heatHT(mid_1+1:summer);
    cost_export_heatHT(summer+1:mid_2)    = weight_mid(1).*cost_export_heatHT(summer+1:mid_2);
    cost_export_heatHT(mid_2+1:winter_2)  = weight_winter(1).*cost_export_heatHT(mid_2+1:winter_2);
    
    cost_export_heatLT(1:winter_1)        = weight_winter(2).*cost_export_heatLT(1:winter_1);
    cost_export_heatLT(winter_1+1:mid_1)  = weight_mid(2).*cost_export_heatLT(winter_1+1:mid_1);
    cost_export_heatLT(mid_1+1:summer)    = weight_summer(2).*cost_export_heatLT(mid_1+1:summer);
    cost_export_heatLT(summer+1:mid_2)    = weight_mid(2).*cost_export_heatLT(summer+1:mid_2);
    cost_export_heatLT(mid_2+1:winter_2)  = weight_winter(2).*cost_export_heatLT(mid_2+1:winter_2);
    
% Parameters for the heat recovery model

    T_in      = 57;                               % Inlet temperature cooling water to HEX    [°C]
    T_out     = 62;                               % Temperature cooling water to applications [°C]
    T_HEX     = 64;                               % Outlet temperature cooling water from HEX [°C]
    c_p       = 4.186;                            % Specific heat capacity cooling water      [kJ/kg/K]
    P_th_max  = (1 - eff_e_nom) * S_e_max;        % Maximum heat recovered                    [kW]
    m_cw_max  = P_th_max/(c_p * (T_out - T_in));  % Cooling water maximum mass flow           [kg/s]
    T_LT_in   = 26;                               % LT DHN water inlet temperature HEX        [°C]
    T_LT_out  = 36;                               % LT DHN water outlet temperature HEX       [°C]
    T_HT_out  = 66;                               % HT DHN water outlet temperature HP        [°C]
    T_log     = ( (T_out - T_LT_out) - (T_in - T_LT_in) ) / ...
       log((T_out - T_LT_out)/(T_in - T_LT_in));  % Logarithmic mean temp difference in HEX
    U_HEX     = 2;                                % Overall heat transfer coeff HEX           [kW/m2/K]
    A_HEX_max = P_th_max / (U_HEX * T_log);       % Maximum heat exchanger surface            [m2]
    COP_carnot= T_HT_out / (T_HT_out - T_out);    % Maximum coefficient of performance HP
    COP       = 0.5 * COP_carnot;                 % Real COP, as in Tiktak


%% Define the optimization problem and the optimization variables

sizingprob = optimproblem;

%% DESIGN VARIABLES

% sizing design variables

    % Electrolyzer size in [W]
    S_e            = optimvar('S_e','LowerBound',S_e_min,'UpperBound',S_e_max);
    % PV area in [m2]
    Area_PV        = optimvar('Area_PV','LowerBound',Area_PV_min,'UpperBound',Area_PV_max);
    % battery capacity in [Wh]
    C_b            = optimvar('C_b','LowerBound',0,'UpperBound',C_b_max);
    % Heat pump size in [W]
    S_HP           = optimvar('S_HP','LowerBound',0,'UpperBound',P_th_max);
    % Heat Exchanger size in m2
    A_HEX          = optimvar('A_HEX','LowerBound',0,'UpperBound', A_HEX_max);
    % cost and sizing parameters for the the compressor
    cost_c         = optimvar('cost_c','LowerBound',0,'UpperBound',cost_c_max);
    y_c1           = optimvar('y_c1','Type','integer','LowerBound',0,'UpperBound',1);
    y_c2           = optimvar('y_c2','Type','integer','LowerBound',0,'UpperBound',1);
    S_cy1          = optimvar('S_cy1','LowerBound',0,'UpperBound',S_c_max);
    S_cy2          = optimvar('S_cy2','LowerBound',0,'UpperBound',S_c_max);
    % cost and sizing parameters for the electrolyzer
    cost_e         = optimvar('cost_e','LowerBound',0,'UpperBound',cost_e_max);
    y_e1           = optimvar('y_e1','Type','integer','LowerBound',0,'UpperBound',1);
    y_e2           = optimvar('y_e2','Type','integer','LowerBound',0,'UpperBound',1);
    S_ey1          = optimvar('S_ey1','LowerBound',0,'UpperBound',S_e_max);
    S_ey2          = optimvar('S_ey2','LowerBound',0,'UpperBound',S_e_max);

% operational design variables
    
    % power consumption for the electrolyzer in [W]
    P_e            = optimvar('P_e',nHours,'LowerBound',0,'UpperBound',S_e_max);
    % Heat generated by the electrolyzer in [W]
    Qdot_H2        = optimvar('Qdot_H2',nHours,'LowerBound',0,'UpperBound',S_e_max); 
    % artificial variables for electrolzyer operation implementation
    P_ElyOn        = optimvar('P_ElyOn',nHours,'LowerBound',0,'UpperBound',S_e_max);
    ElyOn          = optimvar('ElyOn',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    startup        = optimvar('startup',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % Battery energy content
    E_b            = optimvar('E_b',nHours,'LowerBound',0,'UpperBound',C_b_max);                      % realistic boundaries for batteries IRL
    % Battery charging power in [W]
    P_b_ch         = optimvar('P_b_ch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    charging_on    = optimvar('charging_on',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % Battery discharging power in [W]
    P_b_disch      = optimvar('P_b_disch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    discharging_on = optimvar('discharging_on',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % imported power in [W]
    P_imp          = optimvar('P_imp',nHours,'LowerBound',0);
    % exported power in [W]
    P_exp          = optimvar('P_exp',nHours,'LowerBound',0);
    % mass flow rate of cooling water in [kg/s]
    m_cw           = optimvar('m_cw',nHours,'LowerBound',0,'UpperBound',m_cw_max);
    m_cw_HT        = optimvar('m_cw_HT',nHours,'LowerBound',0,'UpperBound',m_cw_max);
    m_cw_LT        = optimvar('m_cw_LT',nHours,'LowerBound',0,'UpperBound',m_cw_max);
    % thermal power transfer to high-tempearture DHN in [W]
    P_th_HT        = optimvar('P_th_HT',nHours,'LowerBound',0,'UpperBound',P_th_max);
    % thermal power transfer to medium-tempearture DHN in [W]
    P_th_LT        = optimvar('P_th_LT',nHours,'LowerBound',0,'UpperBound',P_th_max);


% derived variables

    % PV generated power
    P_PV      = irradiance.*eff_cell*Area_PV/1000;                           % [kW]
    P_PV_peak = 1000*eff_PV*Area_PV/1000;                                    % [kW]
    
    % hydrogen flow calculation
    mdot_H2  = Qdot_H2 * deltat / HHV;          % Mass flow produced hydrogen    [kg/h]
    m_H2_day = sum(reshape(mdot_H2,24,days));   % Mass daily produced hydrogen   [kg/day] sum over the day and vector shaped into a matrix
    P_c = mm_c(3) * mdot_H2;                    % Compressor power               [kW]    
    mdot_H2_nom = eff_e_nom * S_e / HHV;        % Nominal mass flow hydrogen     [kg/s] 
    S_c         = mdot_H2_nom * L_is_c / eff_c; % Size compressor                [kW]
    
    % C-rate of battery
    P_b_lim = C_b / 3600;       

%% CONSTRAINTS

% energy balances

    sizingprob.Constraints.massH2min    = m_H2_day >= H2_demand';
    sizingprob.Constraints.massH2max    = m_H2_day <= 1.1.*H2_demand';
    sizingprob.Constraints.EnBalance    = (P_e + P_th_HT/COP + P_c + P_b_ch/eff_b + P_exp) == (P_PV + eff_b*P_b_disch + P_imp); 
    
    % BATTERY
    sizingprob.Constraints.Batt_charg_fromPV      = P_b_ch/eff_b <= P_PV; % Battery only charged from PV
    sizingprob.Constraints.NoSimultaneousChDisch  = discharging_on + charging_on <= ones(nHours,1);
    sizingprob.Constraints.PowerBatt_ch_0         = P_b_ch(1) == 0; 
    sizingprob.Constraints.PowerBatt_disch_0      = P_b_disch(1) == 0; 
    
    sizingprob.Constraints.Ch_on1        = P_b_ch <= P_b_max * charging_on;
    sizingprob.Constraints.Ch_on2        = P_b_ch <= P_b_lim;
    sizingprob.Constraints.Disch_on1     = P_b_disch <= P_b_max * discharging_on;
    sizingprob.Constraints.Disch_on2     = P_b_disch <= P_b_lim;
    sizingprob.Constraints.E_b           = E_b(idxHr2ToEnd) - E_b(idxHr2ToEnd-1) == P_b_ch(idxHr2ToEnd)*deltat - P_b_disch(idxHr2ToEnd)*deltat;
    sizingprob.Constraints.E_b_cont      = E_b(1) == E_b(end); 
    sizingprob.Constraints.E_b_max       = E_b <= C_b;
    sizingprob.Constraints.E_b_min       = E_b >= 0.2 * C_b;
    sizingprob.Constraints.startupConst  = -ElyOn(idxHr2ToEnd-1) + ElyOn(idxHr2ToEnd) - startup(idxHr2ToEnd) <= 0;
    
    % ELECTROLYSER (EFFICIENCY)
    sizingprob.Constraints.MaxPowerEly   = P_e <= P_ElyOn;
    sizingprob.Constraints.MinPowerEly   = P_e >= 0.2 * P_ElyOn;
    sizingprob.Constraints.PowerEly1     = P_ElyOn <= S_e_max * ElyOn;
    sizingprob.Constraints.PowerEly3     = P_ElyOn <= S_e;
    sizingprob.Constraints.PowerEly4     = P_ElyOn >= S_e - S_e_max.*(ones(nHours,1) - ElyOn);
    
    if N_bp==1
    sizingprob.Constraints.eff_elec      = Qdot_H2 == eff_e_nom * P_e;
    elseif N_bp==2
    sizingprob.Constraints.PWA_Q1        = Qdot_H2 <= mm_elec(1)*P_e + qq_elec(1)*P_ElyOn;
    end

% heat recovery model

    sizingprob.Constraints.PEM_outlet    = m_cw == (P_e-Qdot_H2) / (c_p * (T_HEX - T_in)); % eff_thermal no longer used, now using losses based on van der Roest
    sizingprob.Constraints.cooling       = m_cw >= m_cw_HT + m_cw_LT;
    sizingprob.Constraints.LTheat        = P_th_LT == m_cw_LT * c_p * (T_out - T_in);
    sizingprob.Constraints.HTheat        = P_th_HT == m_cw_HT * c_p * (T_out - T_in);
    sizingprob.Constraints.HEXsize       = P_th_LT / (U_HEX * T_log) <= A_HEX;
    sizingprob.Constraints.HPsize        = P_th_HT <= S_HP;
    sizingprob.Constraints.LTdemand      = P_th_LT <= heat_dem_LT; 
    sizingprob.Constraints.HTdemand      = P_th_HT <= heat_dem_HT;
    

% cost curves

    sizingprob.Constraints.Comp_cy11     = S_cy1 <= S_c_max*y_c1;
    sizingprob.Constraints.Comp_cy13     = S_cy1 <= S_c;
    sizingprob.Constraints.Comp_cy14     = S_cy1 >= S_c - S_c_max*(1-y_c1);
    sizingprob.Constraints.Comp_cy21     = S_cy2 <= S_c_max*y_c2;
    sizingprob.Constraints.Comp_cy23     = S_cy2 <= S_c;
    sizingprob.Constraints.Comp_cy24     = S_cy2 >= S_c - S_c_max*(1-y_c2);
    
    sizingprob.Constraints.Comp_cost     = cost_c == aa_c(1)*S_cy1 + bb_c(1)*y_c1 + aa_c(2)*S_cy2 + bb_c(2)*y_c2;
    sizingprob.Constraints.Comp_cost2    = y_c1 + y_c2 == 1;
    % sizingprob.Constraints.Comp_cost3    = xx_c1*y_c1 + xx_c2*y_c2 <= S_c;
    sizingprob.Constraints.Comp_cost4    = S_c <= xx_c2*y_c1 + xx_c3*y_c2;
    
    sizingprob.Constraints.Elec_ey11     = S_ey1 <= S_e_max*y_e1;
    sizingprob.Constraints.Elec_ey13     = S_ey1 <= S_e;
    sizingprob.Constraints.Elec_ey14     = S_ey1 >= S_e - S_e_max*(1-y_e1);
    sizingprob.Constraints.Elec_ey21     = S_ey2 <= S_e_max*y_e2;
    sizingprob.Constraints.Elec_ey23     = S_ey2 <= S_e;
    sizingprob.Constraints.Elec_ey24     = S_ey2 >= S_e - S_e_max*(1-y_e2);
    
    sizingprob.Constraints.Elec_cost     = cost_e == aa_e(1)*S_ey1 + bb_e(1)*y_e1 + aa_e(2)*S_ey2 + bb_e(2)*y_e2;
    sizingprob.Constraints.Elec_cost2    = y_e1 + y_e2 == 1;
    % sizingprob.Constraints.Elec_cost3    = xx_e1*y_e1 + xx_e2*y_e2 <= S_e;
    sizingprob.Constraints.Elec_cost4    = S_e <= xx_e2*y_e1 + xx_e3*y_e2;

%% OBJECTIVE FUNCTION

cost_inst    = (P_PV_peak * UP_PV * ann_PV + cost_e * ann_e + C_b/3600 * UP_b * ann_b + ...
    S_HP * UP_HP * ann_HP + (A_HEX*UP_HEX + Fixed_HEX) * ann_HEX + cost_c * ann_c + ...
    UP_storage * mass_H2_day_obj * ann_storage + UP_disp * ann_disp + P_refr * UP_refr * ann_refr )/1000; % [kEUR/y]
cost_op      = sum(P_imp.*cost_el)/1000 - sum(P_th_HT.*cost_export_heatHT)/1000 - ...
    sum(P_th_LT.*cost_export_heatLT)/1000 - sum(P_exp.*cost_export_el)/1000; % [kEUR/y]

cost_maint   = (maint_PV * P_PV_peak * UP_PV + maint_e * cost_e + maint_b * C_b/3600 * UP_b * ann_b +...
    maint_HP * S_HP * UP_HP * ann_HP + maint_HEX * (A_HEX*UP_HEX + Fixed_HEX) * ann_HEX + ...
    maint_c * cost_c + maint_disp * UP_disp + maint_refr * P_refr * UP_refr)/1000; % [kEUR/y]

cost_startup = sum(startup*cost_startup_e)/1000;

cost = cost_inst + cost_op + cost_maint + cost_startup;

% set objective
sizingprob.Objective = cost;

%% Solve optimization problem

intcon = [];
options = optimoptions('intlinprog','MaxTime',MaxSimTime);
[solution,fval,reasonSolverStopped] = solve(sizingprob,'Options',options);

%% show problem
% show(sizingprob);

%% Post-processing and results overview

cost_total   = evaluate(cost, solution);
cost_startup = evaluate(cost_startup, solution);
mdot_H2      = evaluate(mdot_H2, solution);
VALCOH         = (cost_total - cost_startup) * 1000 / sum(mdot_H2);

Area_PV_opt   = solution.Area_PV;                                          % [m2]
P_PV_opt      = irradiance.*eff_cell.*solution.Area_PV./1000;              % [kW]
P_e_opt       = solution.P_e;                                              % [W]
P_b_disch_opt = solution.P_b_disch;                                        % [W]
P_b_ch_opt    = solution.P_b_ch;                                           % [W]
P_imp_opt     = solution.P_imp;                                            % [W]
P_exp_opt     = solution.P_exp;                                            % [W]
P_th_opt      = solution.P_e - solution.Qdot_H2;                           % [W]
P_th_av_opt   = solution.m_cw * c_p * (T_out - T_in);                      % [W]
P_th_LT_opt   = solution.P_th_LT;                                          % [W]
P_th_HT_opt   = solution.P_th_HT;                                          % [W]
SOC_opt       = solution.E_b/solution.C_b;                                 % [-]
S_e           = solution.S_e;                                              % [W]
P_PV_peak = evaluate(P_PV_peak, solution);                                 % [W]
S_HP_opt  = solution.S_HP;                                                 % [W]
A_HEX_opt = solution.A_HEX;                                                % [m2]
C_b_opt   = solution.C_b;                                                  % [Wh]

cost_inst = evaluate(cost_inst, solution);                                 % kEUR/y
cost_op   = evaluate(cost_op, solution);                                   % kEUR/y
cost_maint= evaluate(cost_maint, solution);                                % kEUR/y
cost_b    = C_b_opt/3600*UP_b* ann_b /1000;                                 
cost_e    = solution.cost_e * ann_e/1000;
cost_PV   = P_PV_peak*UP_PV* ann_PV/1000;
cost_HP   = S_HP_opt*UP_HP* ann_HP/1000;
cost_HEX  = (A_HEX_opt*UP_HEX + Fixed_HEX) * ann_HEX/1000;
cost_c    = solution.cost_c * ann_c/1000;
cost_th   = - sum(solution.P_th_HT.*cost_export_heatHT)/1000 - sum(solution.P_th_LT.*cost_export_heatLT)/1000;
cost_imp  = sum(solution.P_imp.*cost_el)/1000;
cost_exp  = -sum(solution.P_exp.*cost_export_el)/1000;
Share_PV  = 1-(sum(solution.P_imp)./sum(solution.P_e));
started   = sum(solution.startup);

E_ch  = sum(solution.P_b_ch)/10^3;
E_e   = sum(solution.P_e)/10^3;
E_c   = sum(evaluate(P_c,solution))/10^3; 
E_HP  = sum(solution.P_th_HT./COP)/10^3;
E_exp = sum(solution.P_exp)/10^3;
E_PV    = sum(P_PV_opt)/10^3;
E_disch = sum(solution.P_b_disch)/10^3;
E_imp   = sum(solution.P_imp)/10^3;
E_consumed = E_ch + E_e + E_c + E_HP + E_exp;
E_supplied = E_PV + E_disch + E_imp;
Q_PEM = sum(solution.P_e - solution.Qdot_H2)/10^3;
Q_HP  = sum(solution.P_th_HT)/10^3;
Q_HEX = sum(solution.P_th_LT)/10^3;
Q_amb = Q_PEM - Q_HP - Q_HEX;

Year_char = cellfun(@(v)v(1:4),Input.time(1),'UniformOutput',false);
Year_num  = str2double(Year_char);
OptimalSolution_Summary = [Year_num, days, mass_H2_day_obj, C_b_opt/3600, S_e, P_PV_peak, S_HP_opt, A_HEX_opt, VALCOH, cost_total, cost_inst, cost_maint, cost_op, cost_e, cost_b, cost_HP, cost_HEX, cost_PV, cost_c, cost_th, cost_imp, cost_exp, UP_b, UP_PV, ...
    E_ch, E_e, E_c, E_HP, E_exp, E_PV, E_disch, E_imp, Q_PEM, Q_HP, Q_HEX, Q_amb, Share_PV, MaxSimTime, N_bp, started];

%% Figures - single week plots

% functions directory
addpath(path_functions);

% for a specific week between start and finish
start=7*26*24;
finish=start+24*7+1;

SelectedWeek_SOC(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,S_e,SOC_opt,P_b_ch_opt,P_b_disch_opt)
movegui('center');
 
SelectedWeek_WHR(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,P_th_HT_opt,P_th_LT_opt,solution.ElyOn,P_th_opt,P_th_av_opt)
movegui('east');

SelectedWeek_WHR_ext(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,P_th_HT_opt,P_th_LT_opt,solution.ElyOn)
movegui('west');

