%%      PV CALCULATIONS
% Sources: De Soto W et al. (2006); Dubey et al. (2013); Sun V et al. (2020)

eff_PV_ref  = 0.15;                
T_PV_ref    = 25;
beta_PV_ref = 0.004;
gamma_PV    = 0.12;
NOCT        = 45;
T_cell      = T_amb + (NOCT - 20) * irradiance / 800;

eff_cell    = zeros(length(irradiance),1);

for i = 1:length(irradiance) 
    if irradiance(i) == 0
        eff_cell(i) = 0;

    else
        eff_cell(i)    = eff_PV_ref * ( 1 - beta_PV_ref * (T_cell(i) - T_PV_ref) + ...
            gamma_PV * log10(irradiance(i)));
    end
end


%%      ELECTROLYSER CALCULATIONS
% Using Paolo Gabrielli's script, following is found to describe the
% efficiency of the electrolyser, using the HHV of hydrogen (iso LHV in
% previous script)
% for electrolyser of size 400kW, normalized

if N_bp == 1
    mm_elec = 0.6269;
    qq_elec = 7.8669;
    qq_elec = qq_elec ./ 400; % normalizing the intersection with y-axis, see Paolo Gabrielli

elseif N_bp == 2
    x_bp_val = [77.6555; 428.5104];
    y_bp_val = [56.5484; 276.4963];
    mm_elec = 0.6269;
    qq_elec = 7.8669;
    qq_elec = qq_elec ./ 400; % normalizing the intersection with y-axis, see Paolo Gabrielli
end


%% COMPRESSOR CURVE
% Sources: thesis (Josien de Koning and) Timo Laaksonlaita

mm_c = [4.8756; 6.5348; 2.5470];
med_c = 4.078;                          % Midway H2 flow compressor        [kg/h]
