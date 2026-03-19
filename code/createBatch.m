%% createBatch.m
% Creates a batch structure with a given time duration. 
% A batch is a collection of channels.
% Channels hold the time history for a state/manipulated variable. They
% also contain a description, units of measurement and time vector.
%
% [X] = createBatch(h, T)
%
% h - sample time
% T - batch duration    
function [X] = createBatch(h, T)
    % dummy values
    t = zeros(T/h,1);
    y = zeros(T/h,1);
    
    % pensim manipulated variables
    X.Fg = createChannel('Aeration rate','L/h','h',t,y);
    X.RPM = createChannel('Agitator RPM','RPM','h',t,y);
    X.Fs = createChannel('Sugar feed rate','L/h','h',t,y);
  %  X.Tc = createChannel('Temperature of cooling water','K','h',t,y);
    X.sc = createChannel('Substrate feed concentration','g/L','h',t,y);
    X.abc = createChannel('Acid/base feed concentration','moles','h',t,y);
    X.Fa = createChannel('Acid flow rate','L/h','h',t,y);
    X.Fb = createChannel('Base flow rate','L/h','h',t,y);
    X.Fc = createChannel('Heating/cooling water flow rate','L/h','h',t,y);
    X.Fh = createChannel('Heating water flow rate', 'L/h','h',t,y);

    % indpensim manipulated variables
    X.Fw = createChannel('Water for injection/dilution','L/h','h',t,y);
    X.pressure = createChannel('Air head pressure','bar','h',t,y);
%     X.viscosity = createChannel('Viscosity','centPoise','h',t,y);
    X.Fremoved = createChannel('Dumped broth flow','L/h','h',t,y);
    
    % pensim states
    X.S = createChannel('Substrate concentration','g/L','h',t,y);
    X.DO2 = createChannel('Dissolved oxygen concentration','mg/L','h',t,y);
    X.X = createChannel('Biomass concentration','g/L','h',t,y);
    X.P = createChannel('Penicillin concentration','g/L','h',t,y);
    X.V = createChannel('Vessel Volume','L','h',t,y);
    X.Wt = createChannel('Vessel Weight','Kg','h',t,y);
%     X.CO2 = createChannel('Carbon dioxide concentration','mmole/L','h',t,y);
    X.pH = createChannel('pH','pH','h',t,y);
    X.T = createChannel('Temperature','K','h',t,y);
    X.Q = createChannel('Generated heat','kJ','h',t,y);
  
    % indpensim states
    X.a0 = createChannel('type a0 biomass concentration','g/L','h',t,y);
    X.a1 = createChannel('type a1 biomass concentration','g/L','h',t,y);
    X.a3 = createChannel('type a3 biomass concentration','g/L','h',t,y);
    X.a4 = createChannel('type a4 biomass concentration','g/L','h',t,y);
    X.n0 = createChannel('state n0','-','h',t,y);
    X.n1 = createChannel('state n1','-','h',t,y);
    X.n2 = createChannel('state n2','-','h',t,y);
    X.n3 = createChannel('state n3','-','h',t,y);
    X.n4 = createChannel('state n4','-','h',t,y);
    X.n5 = createChannel('state n5','-','h',t,y);
    X.n6 = createChannel('state n6','-','h',t,y);
    X.n7 = createChannel('state n7','-','h',t,y);
    X.n8 = createChannel('state n8','-','h',t,y);
    X.n9 = createChannel('state n9','-','h',t,y);
    X.nm = createChannel('state nm','-','h',t,y);
    X.phi0 = createChannel('state phi0','-','h',t,y);
    X.CO2outgas = createChannel('carbon dioxide percent in off-gas','%','h',t,y);
    X.Culture_age = createChannel('Cell culture age','h','h',t,y);
    X.Fpaa  =   createChannel('PAA flow','PAA flow (L/h)','h',t,y);
    X.PAA =   createChannel('PAA concentration','PAA (g L^{-1})','h',t,y);
    X.PAA_offline =   createChannel('PAA concentration offline','PAA (g L^{-1})','h',t,y);
    X.Foil =   createChannel('Oil flow','L/hr','h',t,y);
    X.NH3 =   createChannel('NH_3 concentration','NH3 (g L^{-1})','h',t,y);
    X.NH3_offline =   createChannel('NH_3 concentration off-line','NH3 (g L^{-1})','h',t,y);
    X.OUR =   createChannel('OUR','OUR (g min^{-1})','h',t,y);
    X.O2 =   createChannel('Oxygen in percent in off-gas','O2  (%)','h',t,y);
    X.mup =   createChannel('Specific growth rate of Penicillin','mu_P (h^{-1})','h',t,y);
    X.mux =   createChannel('Specific growth rate of Biomass','mu_X (h^{-1})','h',t,y);
    X.P_offline  = createChannel('Offline Penicillin concentration','P(g L^{-1})','h',t,y);
    X.X_CER =  createChannel('Biomass concentration from CER','g min^{-1}','h',t,y);
    X.X_offline = createChannel('Offline Biomass concentratio','X(g L^{-1})','h',t,y);
    X.CER =  createChannel('Carbon evolution rate','g/h','h',t,y);
    X.mu_X_calc =  createChannel('Biomass specific growth rate ','hr^{-1}','h',t,y);
    X.mu_P_calc =  createChannel('Penicillin specific growth rate ','hr^{-1}','h',t,y);
    X.F_discharge_cal = createChannel('Discharge rate','L hr^{-1}','h',t,y); 
    X.NH3_shots = createChannel('Ammonia shots','kgs','h',t,y);
    X.OUR =   createChannel('Oxygen Uptake Rate ','(g min^{-1})','h',t,y);
    X.CO2_d =   createChannel('Dissolved CO_2  ','(mg L^{-1})','h',t,y);
    X.Viscosity =  createChannel('Viscosity ','centPoise','h',t,y);
    X.Viscosity_offline = createChannel('Viscosity','centPoise','h',t,y);
    X.Fault_ref = createChannel('Fault reference', 'Fault ref', 'h', t, y);
    X.Control_ref = createChannel('0 - Recipe driven 1 - Operator controlled', 'Control ref', 'Batch number', t, y);
    X.PAT_ref = createChannel('1- No Raman spec, 1-Raman spec recorded,2-PAT control', 'PAT ref', 'Batch number', t,y);
    X.Batch_ref = createChannel('Batch reference','Batch ref', 'Batch ref', t,y);
  
end
%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037

	