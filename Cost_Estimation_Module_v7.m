clc
clear
cd 'C:\Users\Tirtho Biswas\Documents\MATLAB'
load Model_Results.mat

Dim_op_tab = size(Model_profile_Annual);
Len_op_tab = Dim_op_tab(1,1);
Dummy_var_op_tab = ones(Len_op_tab,1);
%% Exchange Rates

%Exchange Rate Cur to INR[ www.xe.com (Accessed on 24th October 2018)]
Ex_rate_Eur = 84.0211;
Ex_rate_USD = 73.2640;

%% OPEX Cost Variables

%Iron ore pellets price per tonne [https://www.steelmint.com/pellets-prices-indian]
%Iron ore fines can be used directly without agglomeration - [https:https://www.outotec.com/products/direct-and-smelting-reduction/circored-hydrogen-based-reduction/]
%Iron ore pellet price = 7800 INR/tonne
%https://www.business-standard.com/article/economy-policy/iron-ore-price-surge-up-to-17-as-steel-prices-slide-5-118072600991_1.html
%Iron ore fine price = 2400 INR/tonne
%GST rate for iron ore and concentrate 5%
%Ore transport price = USD 7.60 per tonne [BF-BOF costs https://www.steelonthenet.com/cost-bof.html]
ore_price = 2400*1.05 + 7.60*Ex_rate_USD;
%Steel Scrap price per tonne [median taken from ASI 2014-15]
%GST rate 18%
scrap_price = 26379*1.18;
%Flux(limestone) price per tonne [median taken from ASI 2014-15]
%GST 5%
flux_price = 1297*1.05;
%O2 sale price per tonne [https://repository.upenn.edu/cgi/viewcontent.cgi?article=1080&context=cbe_sdr]
%GST rate = 12%
o2_price = 40*73.69*1.12;
%Share of O2 sale
o2_sale_percent = 0.8;
%Grid electricity cost per unit [taken from Guj:http://www.ugvcl.com/petition/Tariff_Schedule.pdf]
%Demand price = INR 475 per month per KVA (in excess of 1000 KVA)
Model_profile_Annual.grid_demand_cost = 475*(Model_profile_Annual{:,{'EAF_power_stoic_MW'}}*1000/(0.9*12))*12;
grid_unit_cost = 4.30;
%Solar LCOE cost [ min 1.5 - max 2.5]
solar_LCOE = 2.3;
%Wind LCOE cost [min 2.5 - max 3.5]
wind_LCOE = 3.4;
%RE price
Model_profile_Annual.re_price = wind_LCOE*Model_profile_Annual.Supply_share_wind + solar_LCOE*Model_profile_Annual.Supply_share_solar;
%re_price = 8;
%Carbon Electrode cost per tonne of electrode(https://www.metalbulletin.com/Article/3824014/INTERVIEW-Breaking-down-the-UHP-graphite-electrode-market.html)
%GST rate 18%
electrode_price = 10000*Ex_rate_USD*1.18;
%Alloy Consumption: 11 kg/tonne of steel
m_alloy = 11*(y_m9*1000)/1000;
%Alloy price
%GST rate 18%
alloy_cost = 1777*Ex_rate_Eur*1.18;
Model_profile_Annual.m_alloy = m_alloy*Dummy_var_op_tab;
%Coke price per tonne [ Assuming USD 250 per tonne]
%GST rate 5% (however for petcoke it is 18%)
cost_coke = 18442.48*1.05;
%Model_profile_Annual.cost_coke = cost_coke*Dummy_var_op_tab;

%% CAPEX cost variables

El_CAPEX_cost = 700*Ex_rate_USD;
%DRI costs per tonne of capacity
DRI_CAPEX_cost = 230*Ex_rate_Eur;
%EAF costs per tonne of capacity
EAF_CAPEX_cost = 184*Ex_rate_Eur;
%H2 storage costs per kg of H2 capacity (USD 14/kwh of stored energy - https://prod.sandia.gov/techlib-noauth/access-control.cgi/2011/114845.pdf)
%GCV (HHV) of H2 = 39.41 kWh/kg
HHV_H2 = 39.41;
%NCV(LHV) of H2 = 33.33 kWh/kg
LHV_H2 = 33.33;
H2store_CAPEX_cost = 14*HHV_H2*Ex_rate_USD;
%SOFC CAPEX costs [USD per KWe](https://www.sainc.com/assets/site_18/files/publications/sa%202015%20manufacturing%20cost%20and%20installed%20price%20of%20stationary%20fuel%20cell%20systems_rev3.pdf)Pg-100 pdf view last page
Fuelcell_CAPEX_cost = 27000*Ex_rate_USD;

%% Plant operational lifetime

%Plant lifetime = 20 years
El_life = 20;
EAF_life = 20;
DRI_life = 20;
Fuelcell_life = 5;
%% Discount Rates

%Refer to discussion from RE team : 80% * 10.5% + 20%*share_wind*16% + 20%*share_solar*14%
%Model_profile_Annual.discount_rate = 0.8*0.105 + 0.2*Model_profile_Annual.Supply_share_wind*0.16 + 0.2*Model_profile_Annual.Supply_share_solar*0.14;
Model_profile_Annual.discount_rate = 0.12*Dummy_var_op_tab;

%% Other variable costs
%Labours costs
%Manhours/tonne of crude steel
labour_man_hours = ((17.4 + 20.1)/2)*(y_m9*1000);
%Hourly cost in India [INR/hr]
labour_hour_cost = 73.69;
labour_cost = labour_man_hours*labour_hour_cost;
Model_profile_Annual.labour_cost = labour_cost*Dummy_var_op_tab;
%By-products [slags, dust, sludges etc. 1 tonne of steel produced in EAF yeilds 200 kg of by-product: https://www.worldsteel.org/en/dam/jcr:1b916a6d-06fd-4e84-b35d-c1d911d18df4/Fact_By-products_2018.pdf]
%byproduct cost = USD 13.07 per tonne
sale_price_byproduct = 13.07*Ex_rate_USD;
%Total by_product generation[tonnes]
m_byproduct = 0.2*y_m9*1000;
revenue_byproduct = sale_price_byproduct*m_byproduct;
revenue_byproduct_pertonne = revenue_byproduct/(y_m9*1000);

%% CAPEX Cost Estimates

Model_profile_Annual.El_CAPEX = El_CAPEX_cost * Model_profile_Annual{:,{'Cap_Elec'}}* 1000;
Model_profile_Annual.DRI_CAPEX = DRI_CAPEX_cost * Model_profile_Annual{:,{'Cap_DRI'}};
Model_profile_Annual.EAF_CAPEX = EAF_CAPEX_cost * Model_profile_Annual{:,{'Cap_EAF'}};
Model_profile_Annual.H2store_CAPEX = H2store_CAPEX_cost* Model_profile_Annual{:,{'m_H2store_stock_annual'}}*24;
%Model_profile_Annual.H2store_CAPEX = H2store_CAPEX_cost* Model_profile_Annual{:,{'m_H2store_stock_annual'}};
Model_profile_Annual.Fuelcell_CAPEX = Fuelcell_CAPEX_cost * Model_profile_Annual{:,{'Max_power_fuelcell_MW'}}*1000;

%Estimating Annual Capital Cost
Model_profile_Annual.El_ACC = Model_profile_Annual.El_CAPEX.*(Model_profile_Annual.discount_rate./(1-((1+Model_profile_Annual.discount_rate).^(-1*El_life))));
Model_profile_Annual.DRI_ACC = Model_profile_Annual.DRI_CAPEX.*(Model_profile_Annual.discount_rate./(1-((1+Model_profile_Annual.discount_rate).^(-1*DRI_life))));
Model_profile_Annual.EAF_ACC = Model_profile_Annual.EAF_CAPEX.*(Model_profile_Annual.discount_rate./(1-((1+Model_profile_Annual.discount_rate).^(-1*EAF_life))));
Model_profile_Annual.H2store_ACC = Model_profile_Annual.H2store_CAPEX.*(Model_profile_Annual.discount_rate./(1-((1+Model_profile_Annual.discount_rate).^(-1*El_life))));
Model_profile_Annual.Fuelcell_ACC = Model_profile_Annual.Fuelcell_CAPEX.*(Model_profile_Annual.discount_rate./(1-((1+Model_profile_Annual.discount_rate).^(-1*Fuelcell_life))));
Model_profile_Annual.Total_ACC = Model_profile_Annual.El_ACC + Model_profile_Annual.DRI_ACC + Model_profile_Annual.EAF_ACC + Model_profile_Annual.H2store_ACC + Model_profile_Annual.Fuelcell_ACC;
%Tot_ACC = (El_ACC + DRI_ACC + EAF_ACC + H2store_ACC)/(y_m9*1000);

%% OPEX Cost Estimates

Model_profile_Annual.resource_cost = (Model_profile_Annual.y_m6/1000)*ore_price + (Model_profile_Annual.y_m8/1000)*scrap_price + (Model_profile_Annual.y_m21/1000)*flux_price + Model_profile_Annual.m_alloy*alloy_cost + Model_profile_Annual.y_mcoke_EAF*cost_coke;
%resource_cost_perTCS = Model_profile_Annual.resource_cost/(y_m9*1000);
Model_profile_Annual.energy_cost = Model_profile_Annual.Energy_El_Annual.*Model_profile_Annual.re_price + Model_profile_Annual.Energy_DRI_Annual.*Model_profile_Annual.re_price + Model_profile_Annual.Energy_EAF_Annual_RE.*Model_profile_Annual.re_price + Model_profile_Annual.Energy_EAF_Annual_Grid*grid_unit_cost + Model_profile_Annual.grid_demand_cost;
Model_profile_Annual.elec_energy_cost = Model_profile_Annual.Energy_El_Annual.*Model_profile_Annual.re_price;
Model_profile_Annual.energy_cost_perTCS = Model_profile_Annual.energy_cost/(y_m9*1000);
%Maintenance and operation costs = 3% of CAPEX
Model_profile_Annual.MO_cost = 0.03 * (Model_profile_Annual.Total_ACC);
Model_profile_Annual.oth_var_cost = (Model_profile_Annual.y_melectrode/1000)*electrode_price + Model_profile_Annual.MO_cost - (Model_profile_Annual.y_m13/1000)*o2_sale_percent*o2_price + Model_profile_Annual.labour_cost;

Model_profile_Annual.OPEX_cost = Model_profile_Annual.resource_cost + Model_profile_Annual.energy_cost + Model_profile_Annual.oth_var_cost;
%Model_profile_Annual.OPEX_cost_perTCS = Model_profile_Annual.resource_cost_perTCS + Model_profile_Annual.energy_cost_perTCS + Model_profile_Annual.oth_var_cost_perTCS;

%% Estimating LCOS

%Cost of steel production in India : http://niti.gov.in/writereaddata/files/document_publication/Need%20for%20a%20new%20Steel%20Policy_NITI%20Website%20Final.pdf
%USD to INR exchange rate : 1 USD = 73.69 INR: 13Oct 2018
avg_cost_TCS_India = (320+340)/2*Ex_rate_USD;
Model_profile_Annual.LCOS_Green_H2_Plant = (Model_profile_Annual.Total_ACC + Model_profile_Annual.OPEX_cost)/(y_m9*1000);
Model_profile_Annual.CUF_Elec = Model_profile_Annual.Energy_El_Annual./(Model_profile_Annual.P1_rated_MW*8760*1000);
save Cost_Results.mat
scatter(Model_profile_Annual.Supply_share_wind, Model_profile_Annual.LCOS_Green_H2_Plant);
writetable(Model_profile_Annual,'Model_profile_Annual.xls');
z = [Model_profile_Annual.Supply_share_wind Model_profile_Annual.CUF_Elec Model_profile_Annual.LCOS_Green_H2_Plant];
surf(z);