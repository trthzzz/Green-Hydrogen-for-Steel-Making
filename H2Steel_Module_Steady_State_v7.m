clc
clear
cd 'C:\Users\Tirtho Biswas\Documents\MATLAB'
%% Description
%Green Hydrogen based DRI-Steel production module
%Process Flow: Electrolysis -- H2 storage -- Condenser/Heat_Exchnager -- Gas_Preheat -- Ore_Preheat -- DRI_Reduction -- HBI_Store/CRDI_Preheat -- Electric_Arc_Furnace
%% Operational Flexibility
%EAF - batch process and tap to tap cycle of 35-40 mins [http://ispatguru.com/understanding-electric-arc-furnace-steel-making-operations/]
%DRI - Instantaneous ramping of 100% feasible [https://www.midrex.cn/assets/user/media/DRI_Products_Brochure_1.29_.15_2.pdf]
%Electrolyser(Elec) - Instantaneous ramping upto 100% has been considered
%% Model Parameters and Key Assumptions
%Input Data - Molar masses [kg]
mm_h2 = 0.002;
mm_o = 0.016;
mm_fe = 0.056;
mm_h2o = 0.018;
mm_fe3o4 = 0.232;
mm_fe2o3 = 0.16;
mm_feo = 0.072;
mm_si02 = 0.060;
mm_caco3 = 0.074;

%Iron to oxyfen ratio for hematite and wuestite
xi_hem = 2/3;
xi_wue = 1;

%module names
%Electrolysis - Elec
%H2 Storage - H2_Store
%Condenser/Heat_Exchanger - Hex
%Gas_Preheat - Gas_Heater
%Ore_Preheat - Ore_Heater
%DRI_Reduction - DRI_Furnace
%HBI_Store/CRDI_Preheat - CRDI_Heater
%Electric_Arc_Furnace - EAF

%Elec
%key Assumptions: No water loss, all water is recovered from shaft off-gas
%Stoichiometric Equations :H2O = H2 + 0.5*O2
ny_h2o_elec = 1;
ny_h2_elec = 1;
ny_o_elec = 1;
%formation enthalpy in J/mol
dHr_el = 242000;
%Efficiency of the Electrolyzer = 72%
eta_el = 0.72;
%Efficiency of the Fuel Cell (SOFC)
eta_fc = 0.60;

%Hex
%Key Assumptions: No heat loss from off-gas between the furnace and exchanger; Exchanger can operate at an efficiency of 75%
eta_hex = 0.75;
%Heat capture by cooling of DRI off gas (H2::H2O) by the exchanger)
%specific heat (gas phase) capacity of Water : Cp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
% For temp range 500 - 1700K
A2_h2o = 30.09200;
B2_h2o = 6.832514;
C2_h2o = 6.793435;
D2_h2o = -2.534480;
E2_h2o = 0.082139;

%For temp range 298 - 500K : liquid phase heat capacity
A1_h2o = -203.6060;
B1_h2o = 1523.290;
C1_h2o = -3196.413;
D1_h2o = 2474.455;
E1_h2o = 3.855326;

%Specific heat (gas phase) capacity of H2: Cp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
%For temp range 298 - 1000K
A1_h2 = 33.066178;
B1_h2 = -11.363417;
C1_h2 = 11.432816;
D1_h2 = -2.772874;
E1_h2 = -0.158558;

%For temp range 1000-2500K
A2_h2 = 18.563083;
B2_h2 = 12.257357;
C2_h2 = -2.859786;
D2_h2 = 0.268238;
E2_h2 = 1.977990;

%DRI
%key Assumptions: No energy loss within the system; DRI is assumed to have a continuous operation. The heat requirement in DRI is separeted into heat requirement to pre-heat ore and gas + heat required for reactions. DRI shaft always running at full capacity
%Operating temperature of DRI = 800C
T_dri = 273.15 + 800;
%Stoichiometric Equations: Fe2O3 + 3*H2 = 2*Fe + 3*H2O; Fe2O3 + H2 = 2*FeO + H2O
ny_fe2o3_dri = 1;
ny_h2_dri = 3;
ny_fe_dri = 2;
ny_h2o_dri = 3;
%formation enthalpy of hematite in J/mol
dHb_fe2o3 = -825500;
%formation enthalpy of water in J/mol
dHb_h2o = -242000;
%reaction enthalpy
dHr_dri = dHb_h2o*ny_h2o_dri - dHb_fe2o3*ny_fe2o3_dri;

ny_fe2o3_dri_wue = 1;
ny_h2_dri_wue = 1;
ny_feo_dri_wue = 2;
ny_h2o_dri_wue = 1;

%DRI-preheating of ore
%Specific heat capacity of Fe(alpha) solid phase is estimated using :Cp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
%Value of A,B,C,D,E are obtained from NIST Chemistry WebBook, SRD 69
%For T =298 to 700K
A1_fe = 18.42868;
B1_fe = 24.64301;
C1_fe = -8.913720;
D1_fe = 9.664706;
E1_fe = -0.012643;

%For T = 700 to 1042K
A2_fe = -57767.65;
B2_fe = 137919.7;
C2_fe = -122773.2;
D2_fe = 38682.42;
E2_fe = 3993.080;

%For T = 1042 to 1100K
A3_fe = -325.8859;
B3_fe = 28.92876;
C3_fe = 0.000000;
D3_fe = 0.000000;
E3_fe = 411.9629;

%DRI-preheating of H2
%Specific heat capacity of H2 (alpha) gaseous phase is estimated using:
%Cp/R_h2 = B + (C-B)*(T/(A+T))^2*[1-((A/(A+T))*(D + E*(T/(A+T)) + F*(T/(A+T))^2 + G*(T/(A+T))^3] - Use same co-efficients from Hex
R_h2 = 4124;

%H2 loading (lambda) of gas feed to the reduction furnace
%H2 loading = moles of H2 present in feed/moles required to reduce unit mass of iron ore in the furnace
%Assuming 50% higher load
lambda_h = 1.5;

%EAF
%key Assumptions: Heat up time of an electric arc furnace is typically very high, EAF always running at full capacity
%Stoichiometric Equations: Fe3O4 = 3*Fe + 4*O; CaCO3 +SiO2 = CaSiO3 + CO2
ny_fe3o4_eaf = 1;
ny_fe_eaf = 3;
ny_o_eaf = 4;
ny_caco3_eaf = 1;
ny_sio3_eaf = 1;
ny_casio3_eaf = 1;
ny_co2_eaf = 1;
%% MODELLING STOIC Operation (Energy and Mass Balances)
%Process Input parameters
y_m9 = input('Enter annual crude steel production in 000 tonnes (must be greater than 400,000) :');
%Converting yearly to daily production capacity (kg per hour)
m9 = (y_m9*10^6)/(365*24);
%y_i8 = input('Enter inert mass fraction in scrap :');
%Coke consumption per tonne of steel consumption[kg]:https://www.researchgate.net/publication/215714915_Study_on_biochar_usage_in_the_electric_arc_furnace/download
mcoke_EAF = 12;
%Inert mass fraction in scrap is assumed to remain at 5%
y_i8 = 0.05;
%y_i6 = input('Enter inert mass fraction in ore pellets :');
%Inert mass fraction in ore is assumed to remain at 5%
y_i6 = 0.06;
alpha = input('Enter degree of metallisation :');
delta = input('Enter scrap blend fraction :');
T_ref = input('Enter the ambient temperature in Celcius :');
T_ref = T_ref + 273.15;
%considering a cross-flow heat exchanger
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EAF - Estimating control variables
%Oxygen content in HBI stream
X_07 = (1-alpha)*mm_o/mm_fe;
%Oxygen content in ore stream
X_06 = 1/xi_hem*mm_o/mm_fe;
%Inert component of HBI stream
X_i6 = y_i6/(1-y_i6);
%Converting oxygen load in ore to mass fraction
y_06 = X_06/(1+ X_06+ X_i6);
%Inert component in scrap stream
X_i8 = y_i8/(1 - y_i8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EAF - Mass Balance
%Iron mass flow to HBI
m7_fe = m9/(1+(delta/(1-delta)*(1-y_i8))*(1+X_07+(1+X_06)*y_i6/(1-y_i6)));
%inert part of m6
m6_i = m7_fe*(1+X_06)*(y_i6/(1-y_i6));
%assuming no iron loss in the DRI reduction
m7 = m7_fe*(1+X_07)+m6_i;
%Scrap flow
m8 = m7*(delta/(1-delta));
%Ore impurities (assuming all impurities to be silicon dioxide)
m_imp = m6_i+m7_fe*(X_07);
moles_imp = m_imp/mm_si02;
m21 = moles_imp*mm_caco3;

%Slag flow (assuming slag = lime + ore impurities + oxygen removed)
m17 = m21+m6_i+m7_fe*(X_07)+m8*y_i8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EAF - Energy Balance
%Starting SEC [J/kg] of EAF is fixed at alpha = 0.94, if alpha and delta are
%varied the SEC is adjusted accordingly
P3_0 = 2.4*10^6;
%Metallisation adjustment factor [J/(kgLS*%metallisation)]
met_adj_fact = 11*3.6*10^3;
%Scrap feed rate adjustment factor [J/(kgLS*%Scrap)]
scrap_adj_fact = 0.85*3.6*10^3;
%EAF power requirement[J/hr]
P3_stoic = (P3_0 + met_adj_fact*(0.94-alpha)*10^2 + scrap_adj_fact*(1-delta)*10^2)*m9;
%EAF power requirement[W]
P3_stoic_W = ((P3_0 + met_adj_fact*(0.94-alpha)*10^2 + scrap_adj_fact*(1-delta)*10^2)*m9)/3600;
%Converting power to yearly units requirments (KWh)
E_y_EAF = (P3_stoic_W*24*365)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRI-HEX - Mass Balance
%Assuming steady state operations and ignoring storage of CRDI
m6_crdi_store = 0;
%Iron content in pellets (assuming 1% loss in iron during reduction in
%furnace)
m6_fe = (m7_fe + 0.01* m7_fe) - m6_crdi_store;
%Estimating ore pellet requirements
m6 = m6_fe*(1+X_06)/(1-y_i6);
m6_stoic = m6;
%Estimating h2 requirement to feed the DRI shaft
m5 = lambda_h*m6_fe*(mm_h2/mm_fe)*(ny_h2_dri/ny_fe_dri);
%Estimating mass of water in shaft off-gas
m10_w = m7_fe*(mm_h2o/mm_fe)*(alpha*ny_h2o_dri/ny_fe_dri+(1-alpha)*(ny_h2o_dri_wue/ny_feo_dri_wue));
%Estimating mass of h2 in shaft off-gas
m10_h = m5-(m7_fe*(mm_h2/mm_fe)*(alpha*(ny_h2_dri/ny_fe_dri) +(1-alpha)*(ny_h2_dri_wue/ny_feo_dri_wue)));
%Total mass of shaft
m10 = m10_w + m10_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Condenser/HEX - Mass Balance
%Assuming no mass loss in the condensing process
m11 = m10_w;
m12 = m10_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modelling HEX energy balance before running DRI energy balance as recovered heat from HEX is used by DRI
%Assuming temperature of off-gas leaving the furnace = 800C
%Assuming temperature of HEX = 70C
T_m12 = 273.15 + 70;
T_m11 = T_m12;

%Estimating heat recovered from cooling of H2 in off-gas by HEX
dt = 1;
T = T_dri;
T_Hex = 273.15 + 70;
h10_h = 0;
while T > T_Hex
    while T > 1000
        % Using gas phase isobaric heat capacity
        Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        h10_h = h10_h + Cp_h*dt;
        %HH10 = m10_h*h10_h;
        T = T-dt;
    end
    %Using gas phase isobaric heat capacity
    Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    h10_h = h10_h + Cp_h*dt;
    %HH10 = HH10+ m10_h*h10_h;
    T = T-dt;
end
HH10_h = h10_h*m10_h;
%Estimating heat recovered from cooling H2O in off-gas by HEX
dt = 1;
T = T_dri;
h10_w = 0;
while T > T_Hex
    while T > 500
        %Using isobaric gas heat phase heat capacity of water
        Cp_w_mol = A2_h2o + B2_h2o*(T/1000) + C2_h2o*((T/1000)^2) + D2_h2o*((T/1000)^3) + E2_h2o/((T/1000)^2);
        Cp_w = Cp_w_mol/mm_h2o;
        h10_w = h10_w + Cp_w*dt;
        T = T-dt;
    end
    %Using isobaric liquid phase heat capacity of water
    Cp_w_mol = A1_h2o + B1_h2o*(T/1000) + C1_h2o*((T/1000)^2) + D1_h2o*((T/1000)^3) + E1_h2o/((T/1000)^2);
    Cp_w = Cp_w_mol/mm_h2o;
    h10_w = h10_w + Cp_w*dt;
    T = T-dt;
end
%Latent heat released while phase change by water [J]
H10_w_latent = m10_w*(2257000/(10^9));
HH10_w = h10_w*m10_w + H10_w_latent;
%HH10_w = h10_w*m10_w;
HH10 = HH10_h + HH10_w;
HH10_hex_stoic = HH10*eta_hex;

%Assuming no heat loss between HEX and the mixing point of H2
T_m3 = T_m12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRI Energy Balance - (Preheating of H2 and Ore + Reduction) Energy Balance[J]
%Preheating of H2
T = T_m3;
dt = 1;
HH5 = 0;
h5 = 0;
while (HH5 < HH10_hex_stoic)&&(T < T_dri)
    while T < 1000
        Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        h5 = h5 + Cp_h*dt;
        T = T + dt;
    end
    Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    h5 = h5 + Cp_h*dt;
    T = T + dt;
end
HH5 = (m5*h5);
%Remaining heat is used for shaft heating[J]
HH_Hex_rec_stoic = HH10_hex_stoic - HH5;
T_m4 = T;

if T_m4 < T_dri
    T = T_m4;
    dt = 1;
    em_5 = 0;
    while T < T_dri
        while T < 1000
        Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        em_5 = em_5 + Cp_h*dt;
        T = T + dt;
        end
    Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    em_5 = em_5 + Cp_h*dt;
    T = T + dt;
    end
   P_m5_stoic = em_5 * m5;
else
    P_m5_stoic = 0;
    if HH_Hex_rec_stoic > 0
         %Recovered heat for shaft heating
         P_Hex_rec_stoic = HH_Hex_rec_stoic;
     end
end

%DRI- Ore preheating from 25C to 800C
T = T_ref;
dt = 1;
h6 = 0;
while T < T_dri
    while T < 700
        Cp_fe_mol = A1_fe + B1_fe*(T/1000) + C1_fe*((T/1000)^2) + D1_fe*((T/1000)^3) + E1_fe/((T/1000)^2);
        Cp_fe = Cp_fe_mol/mm_fe;
        h6 = h6 + Cp_fe*dt;
        T = T + dt;
    end
    while (T > 700 )&&(T < 1042)
        Cp_fe_mol = A2_fe + B2_fe*(T/1000) + C2_fe*((T/1000)^2) + D2_fe*((T/1000)^3) + E2_fe/((T/1000)^2);
        Cp_fe = Cp_fe_mol/mm_fe;
        h6 = h6 + Cp_fe*dt;
        T = T + dt;
    end
    Cp_fe_mol = A3_fe + B3_fe*(T/1000) + C3_fe*((T/1000)^2) + D3_fe*((T/1000)^3) + E3_fe/((T/1000)^2);
    Cp_fe = Cp_fe_mol/mm_fe;
    h6 = h6 + Cp_fe*dt;
    T = T + dt;
end
%Total power required for pre-heating[J/hr]
P_m6_stoic = h6*m6;
%Power required for ore reduction reaction [J/hr]
P_red_stoic = m7_fe*alpha/mm_fe/ny_fe_dri*dHr_dri;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electrolyzer - Mass Balance
%Mixing Point
%Assuming steady state for first run
m_store = 0;
m3 = m5;
m3_stoic = m3;
m2 = m3 - m12 + m_store;
m2_stoic = m2;
m13 = m11/mm_h2o*mm_o;
%Water requirement
m1 = m2 + m13 - m11;

%Electrolyzer - Energy Balance
%Power requirement [J/hr]
P1_stoic = (m2*(1/mm_h2)*dHr_el)/eta_el;
%Energy consumption [KWh]
%E_y_el = P1*24*365/100;
%SEC
EI_el_stoic = P1_stoic/m9;
%Sizing the electrolyzer at steadystate [MW]
Cap_Elec_stoic = P1_stoic/(1000000*3600);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total power requirement[J/hr]
P_total_stoic = P1_stoic + P_red_stoic + P3_stoic - P_Hex_rec_stoic + P_m5_stoic + P_m6_stoic;
%Total power requirement[MW]
P_total_stoic_W = (P1_stoic + P_red_stoic + P3_stoic - P_Hex_rec_stoic + P_m5_stoic + P_m6_stoic)/(3600*10^6);
%Electrolyzer Power[MW]
Elec_power_stoic_MW = P1_stoic/(3600*10^6);
%EAF Power
EAF_power_stoic_MW = P3_stoic/(3600*10^6);
%DRI Power
DRI_power_stoic_MW =(P_red_stoic-P_Hex_rec_stoic+P_m5_stoic + P_m6_stoic)/(3600*10^6);
%Total energy requirement[MWh]
e_y_total_stoic = (P_total_stoic_W*24*365);
%Annual SEC [MWh/tCS]
SEC_stoic = e_y_total_stoic/(y_m9*1000);
save stoic_parameters.mat
%Loading the power profile..refer to the power profile generating module
load RE_profile_gen_module.mat

%% Simulating the hourly production profile for the plant
%%%%%%%%%%%%%%%%%%%%%%%%Generating the optimal RE_profile%%%%%%%%%%%%%%%
REgen_wind_profile_1 = xlsread('RE Dispatch','RE dispatch','D9:D8768');
REgen_wind_profile_2 = xlsread('RE Dispatch','RE dispatch','E9:E8768');
REgen_wind_profile_3 = xlsread('RE Dispatch','RE dispatch','F9:F8768');
REgen_wind_profile_4 = xlsread('RE Dispatch','RE dispatch','G9:G8768');
REgen_wind_profile_installed_capacity_1 = xlsread('RE Dispatch','RE dispatch','D6');
REgen_wind_profile_installed_capacity_2 = xlsread('RE Dispatch','RE dispatch','E6');
REgen_wind_profile_installed_capacity_3 = xlsread('RE Dispatch','RE dispatch','F6');
REgen_wind_profile_installed_capacity_4 = xlsread('RE Dispatch','RE dispatch','G6');
REgen_wind_profile_CUF_1 = (sum(REgen_wind_profile_1)/8760)/REgen_wind_profile_installed_capacity_1;
REgen_wind_profile_CUF_2 = (sum(REgen_wind_profile_2)/8760)/REgen_wind_profile_installed_capacity_2;
REgen_wind_profile_CUF_3 = (sum(REgen_wind_profile_3)/8760)/REgen_wind_profile_installed_capacity_3;
REgen_wind_profile_CUF_4 = (sum(REgen_wind_profile_4)/8760)/REgen_wind_profile_installed_capacity_4;
REgen_wind_profile_CUF = [REgen_wind_profile_CUF_1 REgen_wind_profile_CUF_2 REgen_wind_profile_CUF_3 REgen_wind_profile_CUF_4];
REgen_wind_profile_1MW_1 = REgen_wind_profile_1/REgen_wind_profile_installed_capacity_1;
REgen_wind_profile_1MW_2 = REgen_wind_profile_2/REgen_wind_profile_installed_capacity_2;
REgen_wind_profile_1MW_3 = REgen_wind_profile_3/REgen_wind_profile_installed_capacity_3;
REgen_wind_profile_1MW_4 = REgen_wind_profile_4/REgen_wind_profile_installed_capacity_4;
REgen_wind_profile_1MW = [REgen_wind_profile_1MW_1 REgen_wind_profile_1MW_2 REgen_wind_profile_1MW_3 REgen_wind_profile_1MW_4];
REgen_solar_profile_1 = xlsread('RE Dispatch','RE dispatch','Z9:Z8768');
REgen_solar_profile_2 = xlsread('RE Dispatch','RE dispatch','AA9:AA8768');
REgen_solar_profile_3 = xlsread('RE Dispatch','RE dispatch','AB9:AB8768');
REgen_solar_profile_4 = xlsread('RE Dispatch','RE dispatch','AC9:AC8768');
REgen_solar_profile_installed_capacity_1 = xlsread('RE Dispatch','RE dispatch','Z6');
REgen_solar_profile_installed_capacity_2 = xlsread('RE Dispatch','RE dispatch','AA6');
REgen_solar_profile_installed_capacity_3 = xlsread('RE Dispatch','RE dispatch','AB6');
REgen_solar_profile_installed_capacity_4 = xlsread('RE Dispatch','RE dispatch','AC6');
REgen_solar_profile_CUF_1 = (sum(REgen_solar_profile_1)/8760)/REgen_solar_profile_installed_capacity_1;
REgen_solar_profile_CUF_2 = (sum(REgen_solar_profile_2)/8760)/REgen_solar_profile_installed_capacity_2;
REgen_solar_profile_CUF_3 = (sum(REgen_solar_profile_3)/8760)/REgen_solar_profile_installed_capacity_3;
REgen_solar_profile_CUF_4 = (sum(REgen_solar_profile_4)/8760)/REgen_solar_profile_installed_capacity_4;
REgen_solar_profile_CUF = [REgen_solar_profile_CUF_1 REgen_solar_profile_CUF_2 REgen_solar_profile_CUF_3 REgen_solar_profile_CUF_4];
REgen_solar_profile_1MW_1 = REgen_solar_profile_1/REgen_solar_profile_installed_capacity_1;
REgen_solar_profile_1MW_2 = REgen_solar_profile_2/REgen_solar_profile_installed_capacity_2;
REgen_solar_profile_1MW_3 = REgen_solar_profile_3/REgen_solar_profile_installed_capacity_3;
REgen_solar_profile_1MW_4 = REgen_solar_profile_4/REgen_solar_profile_installed_capacity_4;
REgen_solar_profile_1MW = [REgen_solar_profile_1MW_1 REgen_solar_profile_1MW_2 REgen_solar_profile_1MW_3 REgen_solar_profile_1MW_4];
REgen_installed_profile = 0;
k = 1;
j = 1;
i = 0;
Model_profile_Annual = table;
while i <= 1    
%wind_installed_capacity = ((DRI_power_stoic_MW + Elec_power_stoic_MW)*(1-i))/REgen_wind_profile_CUF(j);
%solar_installed_capacity = ((DRI_power_stoic_MW + Elec_power_stoic_MW)*i)/REgen_solar_profile_CUF(k);
%%wind_installed_capacity = 1310*(1-i);
%%solar_installed_capacity = 368.7*i;
%wind_installed_profile = REgen_wind_profile_1MW(:,j)*wind_installed_capacity;
%solar_installed_profile = REgen_solar_profile_1MW(:,k)*solar_installed_capacity;
%REgen_installed_profile = solar_installed_profile + wind_installed_profile;
Profile_parameters = [ j 1-i k i];
Supply_share_wind = 1-i;
Supply_share_solar = i;
%%%%%%%Running the plant operations from here%%%%%%%%%%%
P1_rated = 0;
m_H2store_stock_recharge = -0.000000000001;
m_H2store_stock_withdraw = 0;
%X = 1;
X = 1;
count = 0;
while m_H2store_stock_recharge < m_H2store_stock_withdraw
     %Sizing Electrolyser capacity
     P1_rated = X*P1_stoic;
     P1_rated_MW = P1_rated/(3600*10^6);
     wind_installed_capacity = ((DRI_power_stoic_MW + P1_rated_MW)*(1-i))/REgen_wind_profile_CUF(j);
     solar_installed_capacity = ((DRI_power_stoic_MW + P1_rated_MW)*i)/REgen_solar_profile_CUF(k);
     wind_installed_profile = REgen_wind_profile_1MW(:,j)*wind_installed_capacity;
     solar_installed_profile = REgen_solar_profile_1MW(:,k)*solar_installed_capacity;
     REgen_installed_profile = solar_installed_profile + wind_installed_profile;
     P1_min = 0;
     Hours = 1;
     Power_profile_elec_annual = zeros(1,2);
     Power_profile_fuelcell_annual = 0;
     m_H2store_stock_wihdraw_profile_annual = 0;
     Op_Hours_El = 0;
     Energy_El_Annual = 0;
     Energy_DRI_Annual = 0;
     Energy_EAF_Annual_Grid = 0;
     Energy_EAF_Annual_RE = 0;
     y_m5 = 0;
     y_m6 = 0;
     y_m13 = 0;
     y_m21 = 0;
     y_m8 = 0;
     energy_curtail_MWh = 0;
     y_melectrode = 0;

while Hours <= 8760
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EAF - Estimating control variables
%Oxygen content in HBI stream
X_07 = (1-alpha)*mm_o/mm_fe;
%Oxygen content in ore stream
X_06 = 1/xi_hem*mm_o/mm_fe;
%Inert component of HBI stream
X_i6 = y_i6/(1-y_i6);
%Converting oxygen load in ore to mass fraction
y_06 = X_06/(1+ X_06+ X_i6);
%Inert component in scrap stream
X_i8 = y_i8/(1 - y_i8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EAF - Mass Balance
%Iron content in mass flow from HBI per hour
m7_fe = m9/(1+(delta/(1-delta)*(1-y_i8))*(1+X_07+(1+X_06)*y_i6/(1-y_i6)));
%inert part of m6
m6_i = m7_fe*(1+X_06)*(y_i6/(1-y_i6));
%assuming no iron loss in the DRI reduction
m7 = m7_fe*(1+X_07)+m6_i;
%Scrap flow
m8 = m7*(delta/(1-delta));
y_m8 = y_m8 + m8;
%Ore impurities (assuming all impurities to be silicon dioxide)
m_imp = m6_i+m7_fe*(X_07);
moles_imp = m_imp/mm_si02;
%Assuming lime feed of 50Kg/tonnes of liquid steel[kg]
m21 = 0.050*m9;
y_m21 = y_m21 + m21;
%Electrode consumption[kg]
melectrode = 0.002*m9;
y_melectrode = y_melectrode + melectrode;
%Slag flow (assuming slag = lime + ore impurities + oxygen removed)
m17 = m21+m6_i+m7_fe*(X_07)+m8*y_i8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EAF - Energy Balance
%Starting SEC [J/kg] of EAF is fixed at alpha = 0.94, if alpha and delta are varied the SEC is adjusted accordingly
P3_0 = 2.4*10^6;
%Metallisation adjustment factor [J/(kgLS*%metallisation)]
met_adj_fact = 11*3.6*10^3;
%Scrap feed rate adjustment factor [J/(kgLS*%Scrap)]
scrap_adj_fact = 0.85*3.6*10^3;
%EAF power requirement[J/hr]
P3 = (P3_0 + met_adj_fact*(0.94-alpha)*10^2 + scrap_adj_fact*(1-delta)*10^2)*m9;
%EAF energy requirement[KWh]
Energy_EAF = P3/(3600*1000);
%Diverting additional power (after DRI and Elec) to EAF [J/hr];
P_Avail = REgen_installed_profile(Hours,:)*(10^6)*(3600);
P1_Avail = P_Avail - (DRI_power_stoic_MW)*(10^6)*(3600);
if P1_Avail > P1_rated
    P3_RE = P1_Avail - P1_rated;
else
    P3_RE = 0;
end
%Energy supply from RE
if P3_RE >= P3
    Energy_EAF_RE = Energy_EAF;
else
    Energy_EAF_RE = P3_RE/(3600*1000);
end
%Energy supply from Grid
Energy_EAF_Grid = Energy_EAF - Energy_EAF_RE;
%Annual Energy requirement[KWh]
Energy_EAF_Annual_RE = Energy_EAF_Annual_RE + Energy_EAF_RE;
Energy_EAF_Annual_Grid = Energy_EAF_Annual_Grid + Energy_EAF_Grid;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRI-HEX - Mass Balance
%Assuming steady state operations and ignoring storage of CRDI
m6_crdi_store = 0;
%Iron content in ore (fines) (assuming 1% loss in iron during reduction in furnace)
m6_fe = (m7_fe + 0.01* m7_fe) - m6_crdi_store;
%Estimating ore ore requirements
m6 = m6_fe*(1+X_06)/(1-y_i6);
y_m6 = y_m6 + m6;
%Estimating h2 requirement to feed the DRI shaft
m5 = lambda_h*m6_fe*(mm_h2/mm_fe)*(ny_h2_dri/ny_fe_dri);
y_m5 = y_m5 + m5;
%Estimating mass of water in shaft off-gas
m10_w = m7_fe*(mm_h2o/mm_fe)*(alpha*ny_h2o_dri/ny_fe_dri+(1-alpha)*(ny_h2o_dri_wue/ny_feo_dri_wue));
%Estimating mass of h2 in shaft off-gas
m10_h = m5-(m7_fe*(mm_h2/mm_fe)*(alpha*(ny_h2_dri/ny_fe_dri) +(1-alpha)*(ny_h2_dri_wue/ny_feo_dri_wue)));
%Total mass of shaft
m10 = m10_w + m10_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Condenser/HEX - Mass Balance
%Assuming no mass loss in the condensing process
m11 = m10_w;
m12 = m10_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modelling HEX energy balance before running DRI energy balance as recovered heat from HEX is used by DRI
%Assuming temperature of off-gas leaving the furnace = 800C
%Assuming temperature of HEX = 70C
T_m12 = 273.15 + 70;
T_m11 = T_m12;

%Estimating heat recovered from cooling of H2 in off-gas by HEX
dt = 1;
T = T_dri;
T_Hex = 273.15 + 70;
h10_h = 0;
while T > T_Hex
    while T > 1000
        % Using gas phase isobaric heat capacity
        Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        h10_h = h10_h + Cp_h*dt;
        %HH10 = m10_h*h10_h;
        T = T-dt;
    end
    %Using gas phase isobaric heat capacity
    Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    h10_h = h10_h + Cp_h*dt;
    %HH10 = HH10+ m10_h*h10_h;
    T = T-dt;
end
HH10_h = h10_h*m10_h;
%Estimating heat recovered from cooling H2O in off-gas by HEX
dt = 1;
T = T_dri;
h10_w = 0;
while T > T_Hex
    while T > 500
        %Using isobaric gas heat phase heat capacity of water
        Cp_w_mol = A2_h2o + B2_h2o*(T/1000) + C2_h2o*((T/1000)^2) + D2_h2o*((T/1000)^3) + E2_h2o/((T/1000)^2);
        Cp_w = Cp_w_mol/mm_h2o;
        h10_w = h10_w + Cp_w*dt;
        T = T-dt;
    end
    %Using isobaric liquid phase heat capacity of water
    Cp_w_mol = A1_h2o + B1_h2o*(T/1000) + C1_h2o*((T/1000)^2) + D1_h2o*((T/1000)^3) + E1_h2o/((T/1000)^2);
    Cp_w = Cp_w_mol/mm_h2o;
    h10_w = h10_w + Cp_w*dt;
    T = T-dt;
end
%Latent heat released while phase change by water [J]
H10_w_latent = m10_w*(2257000/(10^9));
HH10_w = h10_w*m10_w + H10_w_latent;
%HH10_w = h10_w*m10_w;
HH10 = HH10_h + HH10_w;
HH10_hex = HH10*eta_hex;

%Assuming no heat loss between HEX and the mixing point of H2
T_m3 = T_m12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRI - (Preheating of H2 and Ore + Reduction) Energy Balance[J]
%Preheating of H2
T = T_m3;
dt = 1;
HH5 = 0;
h5 = 0;
while (HH5 < HH10_hex)&&(T < T_dri)
    while T < 1000
        Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        h5 = h5 + Cp_h*dt;
        T = T + dt;
    end
    Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    h5 = h5 + Cp_h*dt;
    T = T + dt;
end
HH5 = -1* (m5*h5);
%Remaining heat is used for shaft heating[J]
HH_Hex_rec = HH10_hex - HH5;
T_m4 = T;

if T_m4 < T_dri
    T = T_m4;
    dt = 1;
    em_5 = 0;
    while T < T_dri
        while T < 1000
        Cp_h_mol = A1_h2 + B1_h2*(T/1000) + C1_h2*((T/1000)^2) + D1_h2*((T/1000)^3) + E1_h2/((T/1000)^2);
        Cp_h = Cp_h_mol/mm_h2;
        em_5 = em_5 + Cp_h*dt;
        T = T + dt;
        end
    Cp_h_mol = A2_h2 + B2_h2*(T/1000) + C2_h2*((T/1000)^2) + D2_h2*((T/1000)^3) + E2_h2/((T/1000)^2);
    Cp_h = Cp_h_mol/mm_h2;
    em_5 = em_5 + Cp_h*dt;
    T = T + dt;
    end
   P_m5 = em_5 * m5;
else
    P_m5 = 0;
    if HH_Hex_rec > 0
         %Recovered heat for shaft heating
         P_Hex_rec = HH_Hex_rec;
     end
end

%DRI- Ore preheating from 25C to 800C
T = T_ref;
dt = 1;
h6 = 0;
while T < T_dri
    while T < 700
        Cp_fe_mol = A1_fe + B1_fe*(T/1000) + C1_fe*((T/1000)^2) + D1_fe*((T/1000)^3) + E1_fe/((T/1000)^2);
        Cp_fe = Cp_fe_mol/mm_fe;
        h6 = h6 + Cp_fe*dt;
        T = T + dt;
    end
    while (T > 700 )&&(T < 1042)
        Cp_fe_mol = A2_fe + B2_fe*(T/1000) + C2_fe*((T/1000)^2) + D2_fe*((T/1000)^3) + E2_fe/((T/1000)^2);
        Cp_fe = Cp_fe_mol/mm_fe;
        h6 = h6 + Cp_fe*dt;
        T = T + dt;
    end
    Cp_fe_mol = A3_fe + B3_fe*(T/1000) + C3_fe*((T/1000)^2) + D3_fe*((T/1000)^3) + E3_fe/((T/1000)^2);
    Cp_fe = Cp_fe_mol/mm_fe;
    h6 = h6 + Cp_fe*dt;
    T = T + dt;
end
%Total power required for pre-heating[J]
P_m6 = h6*m6;
%Power required for ore reduction reaction [J]
P_red = m7_fe*alpha/mm_fe/ny_fe_dri*dHr_dri;
%Annual energy consumption by the DRI Furnace [KWh]
Energy_DRI_Annual = Energy_DRI_Annual + (P_red_stoic-P_Hex_rec_stoic+P_m5_stoic + P_m6_stoic)/(3600*10^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electrolyzer - Mass Balance
%Mixing Point
m3 = m5;
%m2 = m3 - m12 + m_H2store_withdraw;
m13 = m11/mm_h2o*mm_o;
y_m13 = y_m13 + m13;

%Mass balance at electrolyser
%m1 = m2 + m13 - m11;

%Electrolyzer - Energy Balance
%Power requirement [J]
%P1_Avail = Available power to the electrolyser from RE_profile[J/hr] after discounting power from DRI and EAF operations;
if P1_Avail >= P1_stoic
    %Electrolyser would produce more hydrogen than steady state and the excess hydrogen would be used to store
    if P1_Avail > P1_rated
        %Curtailing excess power (power higher than the rated capacity)
        P1_curtail = P1_Avail - 0.95*P1_rated;
        P1_Avail = 0.95*P1_rated;
        m2 = (P1_Avail * eta_el * mm_h2 )/dHr_el;
        m_H2store_stock_recharge = m_H2store_stock_recharge + (m12 + m2)-m3;
    else
        m2 = (P1_Avail * eta_el * mm_h2 )/dHr_el;
        m_H2store_stock_recharge = m_H2store_stock_recharge + (m12 + m2)-m3;
    end
    %Since curtailment is done every hour every conversion of J/hr to MW would represent MWh
    energy_curtail_MWh = energy_curtail_MWh + P1_curtail /(3600*10^6);
    %Electrolyser would produce less hydrogen and remaining hydrogen would be withdrawn from storage
else
    if P1_Avail > P1_min
        m2 = (P1_Avail * eta_el * mm_h2 )/dHr_el;
        m_H2store_stock_withdraw = m_H2store_stock_withdraw + m3 - (m12 + m2);
        m_H2store_stock_wihdraw_profile = m3 - (m12 + m2);
        m_H2store_stock_wihdraw_profile_annual = vertcat(m_H2store_stock_wihdraw_profile_annual, m_H2store_stock_wihdraw_profile);
    else
        P1_Avail = 0;
        m2 = (P1_Avail * eta_el * mm_h2 )/dHr_el;
        %If available power falls below the DRI power requirement, a fuel cell is run to meet the demand
        P_fuelcell = DRI_power_stoic_MW*10^6*3600 - P_Avail;
        P_fuelcell_MW = P_fuelcell/(3600*10^6);
        Power_profile_fuelcell_annual = vertcat(Power_profile_fuelcell_annual,P_fuelcell_MW);
        m_H2fuelcell = (P_fuelcell*mm_h2 )/(dHr_el*eta_fc);
        m_H2store_stock_withdraw = m_H2store_stock_withdraw + m3 - (m12 + m2) + m_H2fuelcell;
        m_H2store_stock_wihdraw_profile =  m3 - (m12 + m2) + m_H2fuelcell;
        m_H2store_stock_wihdraw_profile_annual = vertcat(m_H2store_stock_wihdraw_profile_annual, m_H2store_stock_wihdraw_profile);
    end
end
%No of operating hours of the electrolyser
if P1_Avail>0
     Op_Hours_El =  Op_Hours_El +1;
end
%Annual energy consumption by electrolyser in KWh
Energy_El_Annual = Energy_El_Annual + P1_Avail/(3600*1000);
%Creating the Electrolyser energy consumption profile
Power_profile_elec = [Hours P1_Avail];
Power_profile_elec_annual = vertcat(Power_profile_elec_annual, Power_profile_elec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total power requirement[J/hr]
P_total = P1_Avail + P_red + P3 - P_Hex_rec + P_m5 + P_m6;
Hours = Hours + 1;
end
m_H2store_stock_annual = max(m_H2store_stock_wihdraw_profile_annual);
X = X+0.05
count = count + 1
if count == 100
    flag = 1;
    break
end
end
if flag == 1
    break
end
%Maximum power supplied by fuel cell
Max_power_fuelcell_MW = max(Power_profile_fuelcell_annual);
%Electrolyser capacity in MW
Cap_Elec = P1_rated/(3600*10^6);
%Cap in tonnes
Cap_EAF = y_m9*1000;
%Cap in tonnes
Cap_DRI = y_m6/1000;
%Total coke consumption[tonnes]
y_mcoke_EAF = (mcoke_EAF*(y_m9*1000))/1000;
Model_profile = table(Supply_share_wind, Supply_share_solar, EAF_power_stoic_MW, DRI_power_stoic_MW, Elec_power_stoic_MW, P1_rated_MW, P1_rated, m_H2store_stock_annual, y_m6, y_m5, y_m13, y_m21, y_m8, y_melectrode, Energy_EAF_Annual_RE, Energy_EAF_Annual_Grid,Max_power_fuelcell_MW, Energy_El_Annual, Energy_DRI_Annual, y_mcoke_EAF, Cap_Elec, Cap_EAF, Cap_DRI, energy_curtail_MWh); 
Model_profile_Annual = vertcat(Model_profile_Annual, Model_profile);
i = i + 0.1
end
save Model_Results.mat
