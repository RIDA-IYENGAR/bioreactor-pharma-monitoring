%% fctrl_indpensim.m
% Provides manipulated variables to control IndPenSim
%
% X - batch data. DO NOT CHANGE VALUE IN CODE!!!
%     Check pensim3.m for more information on it's contents
% Xd - industrial data
% k - simulation sample
% h - simulation sample period
% T - simulation time duration
%
% u - structure with the values of the manipulated variables
%   u.Fg            aeration rate [L/h]
%   u.Pw            aeration power [W]
%   u.F             substrate feed rate [L/h]
%   u.Tf            substrate feed temperature [K]
%   u.Fa            acid flow rate [L/h]
%   u.Fb            base flow rate [L/h]
%   u.Fc            heating/cooling water flow rate [L/h]
%   u.Fw         water for injection/dilution [L/h]
%   u.pressure      air head pressure [bar] ?? why is this an input?
%   u.viscosity     broth viscosity input [cPoise]
%   u.Fremoved      dumped broth flow [L/h]

function [u X] = fctrl_indpensim(X, Xd, k, h, T, Ctrl_flags)

%% pH controller
% Adding in Temperature sensor fault
pH_sensor_error = 0; % no error
if Ctrl_flags.Faults== 8; 
pH_sensor_error = 0.1; 
Ramp_function = [ 0 0; 200 0;800 pH_sensor_error; 1750 pH_sensor_error];
tInterp = 1:1:1750; 
Ramp_function_interp = interp1(Ramp_function(:,1),Ramp_function(:,2), tInterp,'linear','extrap');
pH_sensor_error = Ramp_function_interp(k);  
u.Fault_ref = 1; % Fault occured
end 
% builds the error history. Samples 1 and 2 are calculated separately 
% because there is only instance of the error available
pH_sp = Ctrl_flags.pH_sp; % pH set-point (-)
if(k==1 || k==2)
    ph_err  = pH_sp - (-log(X.pH.y(1))/log(10))+pH_sensor_error;
    ph_err1  = pH_sp - (-log(X.pH.y(1))/log(10))+pH_sensor_error;
else
    ph_err  = pH_sp - (-log(X.pH.y(k-1))/log(10))+pH_sensor_error;
    ph_err1  =  - (-log(X.pH.y(k-2))/log(10))+pH_sensor_error;
end

% builds the pH history of the current and previous two samples
if(k==1 || k==2)
    ph = -log(X.pH.y(1))/log(10);
    ph1 = -log(X.pH.y(1))/log(10);
    ph2 = -log(X.pH.y(1))/log(10);
elseif (k==3)
    ph = -log(X.pH.y(2))/log(10);
    ph1 = -log(X.pH.y(1))/log(10);
    ph2 = -log(X.pH.y(1))/log(10);
else
    ph = -log(X.pH.y(k-1))/log(10);
    ph1 = -log(X.pH.y(k-2))/log(10);
    ph2 = -log(X.pH.y(k-3))/log(10);
end

% pH has decreased under 0.05 from set-point, add some base solution
if ph_err>=-0.05 
    ph_on_off=1; % former d1 in initial simulation
    % maximum base addition rate is 0.1L/h
 if(k==1)                                                       % 0.005
        Fb = PIDSimple3(X.Fb.y(1), ph_err, ph_err1, ph, ph1, ph2, 0, 225, 8e-2, 4.0000e-05, 8,h);
    else
             Fb = PIDSimple3(X.Fb.y(k-1), ph_err, ph_err1, ph, ph1, ph2, 0, 225, 8e-2, 4.0000e-05, 8,h);     
    end
    Fa = 0;% acid flow rate    
% pH has increased to more than 0.05 from set-point, add some acid. There is a threshold
% level to prevent frequent acid adding
elseif ph_err<=-0.05
      ph_on_off=1;
  if(k==1)
          Fa = PIDSimple3(X.Fa.y(1), ph_err, ph_err1, ph, ph1, ph2, 0, 225, 8e-2, 12.5, 0.125,h);
          Fb=0;
         else
         Fa = PIDSimple3(X.Fa.y(k-1), ph_err, ph_err1, ph, ph1, ph2, 0, 225, 8e-2, 12.5, 0.125,h);
        
        Fb = X.Fb.y(k-1)*0.5;
         end
% pH is on set-point, do nothing
else
    ph_on_off = 0;
     Fb = 0;
     Fa = 0;
end
%% Temperature controller
% Adding in Temperature sensor fault
T_sensor_error = 0; % no error
if Ctrl_flags.Faults== 7; 
T_sensor_error = 0.4; 
Ramp_function = [ 0 0; 200 0;800 T_sensor_error; 1750 T_sensor_error];
tInterp = 1:1: 1750; 
Ramp_function_interp = interp1(Ramp_function(:,1),Ramp_function(:,2), tInterp,'linear','extrap');
T_sensor_error = Ramp_function_interp(k);   
u.Fault_ref = 1; % Fault occured
end 
% builds the error history.  Samples 1 and 2 are calculated separately 
% because there is only instance of the error available
T_sp =  Ctrl_flags.T_sp; 
if(k==1 || k==2)
    temp_err  = T_sp-X.T.y(1)+T_sensor_error;
    temp_err1 = T_sp-X.T.y(1)+T_sensor_error;
else
    temp_err  = T_sp-X.T.y(k-1)+T_sensor_error;
    temp_err1 = T_sp-X.T.y(k-2)+T_sensor_error;
end

% builds the temperature history of current and previous two samples.
if(k==1 || k==2)
    temp = X.T.y(1);
    temp1 = X.T.y(1);
    temp2 = X.T.y(1);
elseif (k==3)
    temp = X.T.y(2);
    temp1 = X.T.y(1);
    temp2 = X.T.y(1);
else
    temp = X.T.y(k-1);
    temp1 = X.T.y(k-2);
    temp2 = X.T.y(k-3);
end

% % % Treshold for heating. Heating is activated only if the temperature drop
% % is more than 1 degree celsius
if temp_err<=0.05
     temp_on_off = 0; % cooling 
    if(k==1)
        [Fc ] = PIDSimple3(X.Fc.y(1), temp_err, temp_err1, temp, temp1, temp2,  0, 1.5e3, -300, 1.6, 0.005, h);
         Fh = 0;
    else
        [Fc] = PIDSimple3(X.Fc.y(k-1), temp_err, temp_err1, temp, temp1, temp2, 0, 1.5e3, -300, 1.6 ,0.005, h);
      Fh = X.Fh.y(k-1)*0.1;
    end
    
    else
     temp_on_off = 1; % heating
    % !gain is negative
    if(k==1)
        [Fh]= PIDSimple3(X.Fc.y(1), temp_err, temp_err1, temp, temp1, temp2, 0, 1.5e3, 50, 0.050, 1, h);
    Fc = 0;
    else
        [Fh] = PIDSimple3(X.Fc.y(k-1), temp_err, temp_err1, temp, temp1, temp2, 0, 1.5e3, 50, 0.050, 1, h);
        Fc = 0; 
         Fc = X.Fc.y(k-1)*0.3;
    end
end 
% necessary for numerical stability
if(Fc < 1e-4)
    Fc = 1e-4;
end
if(Fh<1e-4)
 Fh = 1e-4;
end

%% Sequential Batch control strategy
% If Sequential Batch Control (SBC) = 1, operator controlled
if Ctrl_flags.SBC ==1
Foil =    Xd.Foil.y(k); 
F_discharge  = Xd.F_discharge_cal.y(k);
pressure = Xd.pressure.y(k);
Fpaa = Xd.Fpaa.y(k);
Fw = Xd.Fw.y(k);
viscosity = Xd.viscosity.y(k);
% User input aeration rates
Fg = Xd.Fg.y(k);
Fs = Xd.Fs.y(k);
end
% If Sequential Batch Control (SBC) = 0, controlled using SBC as described
% below, these sequences can be modified for better control
%% SBC - Fs
if Ctrl_flags.SBC ==0 
viscosity = 4;
Recipe_Fs =    [15 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 800 1750 ]; % Sequential batch control standard 
%Recipe_Fs_sp = [8  15 30 75  150 30  37  43  47  51  57   61 65   72 76  80  84  90  136 130 120 ];
Recipe_Fs_sp = [8  15 30 75  150 30  37  43  47  51  57   61 65   72 76  80  84  90  116 90 80 ]; % Adjusted
% Example of modified Fs control strategy
%   Recipe_Fs =    [15 60 80 100 140 160 180 200 220 240 260 280 300 320 340 360 380 400 800 1750 ]; % Sequential batch control standard 
%   Recipe_Fs_sp = [8  15 30 75  300  37  43  47  51  57   61 65   72 76  80  84  90  156 156 120 ];
  for SQ = 1:1:length(Recipe_Fs) % Fs Sequence number 
if k<=Recipe_Fs(SQ)
    Fs =  Recipe_Fs_sp(SQ);
    break
else 
    Fs =  Recipe_Fs_sp(end);
end 
  end
% Add PRBS to substrate flow rate  (Fs)
if Ctrl_flags.PRBS== 1;
if k > 500 && rem(k,100)==0
    rng('shuffle');
random_number = randi(3,1,1);
    noise_factor =15; 
if random_number == 1; 
    random_noise = 0;
elseif random_number == 2;
    random_noise = noise_factor;
else 
    random_noise = -noise_factor;
end 
X.PRBS_noise_addition(k) =random_noise;
end 
if k >475
     Fs = X.Fs.y(k-1);
end
if  k>500 && rem(k,100)==0
Fs = X.Fs.y(k-1) + X.PRBS_noise_addition(end);
end 
else 
    X.PRBS_noise_addition(k) =0;
    
end 
  
%% SBC - Foil
Recipe_Foil =    [20 80 280 300 320 340 360 380 400  1750 ]; % Sequential batch control 
Recipe_Foil_sp = [22 30 35  34  33  32  31  30  29  23 ];
for SQ = 1:1:length(Recipe_Foil) % Foil Sequence number 
% User input aeration rates
if k<=Recipe_Foil(SQ)
    Foil =  Recipe_Foil_sp(SQ);
    break
else 
    Foil =  Recipe_Foil_sp(end);
end 
end 
%% SBC - Fg
Recipe_Fg =    [40 100  200 450 1000 1250  1750 ]; % Sequential batch control 
Recipe_Fg_sp = [30  42  55   60  75   65   60 ];
for SQ = 1:1:length(Recipe_Fg) % Fg Sequence number 
% User input aeration rates
if k<=Recipe_Fg(SQ)
    Fg =  Recipe_Fg_sp(SQ);
    break
else 
    Fg =  Recipe_Fg_sp(end);
end 
end
%% SBC - pressure
Recipe_pres =    [62.5 125  150 200 500 750 1000  1750 ]; % Sequential batch control 
Recipe_pres_sp = [0.6  0.7  0.8 0.9 1.1  1  0.9  0.9 ];
for SQ = 1:1:length(Recipe_pres) % Fg Sequence number 
% User input aeration rates
if k<=Recipe_pres(SQ)
    pressure =  Recipe_pres_sp(SQ);
    break
else 
    pressure =  Recipe_pres_sp(end);
end 
end
% SBC - F_discharge
Recipe_discharge =    [500 510  650 660  750 760  850 860 950 960 1050 1060 1150 1160 1250 1260 1350 1360 1750 ]; % Sequential batch control 
Recipe_discharge_sp = [ 0  4000  0  4000  0  4000 0  4000 0  4000 0  4000 0 4000 0    4000 0    4000  0     0 ];
for SQ = 1:1:length(Recipe_discharge) % Fg Sequence number 
% User input aeration rates
if k<=Recipe_discharge(SQ)
    F_discharge =  -Recipe_discharge_sp(SQ);
    break
else 
    F_discharge =  Recipe_discharge_sp(end);
end 
end
%% SBC - Fw
Recipe_water =    [250 375  750 800 850  1000 1250 1350  1750 ]; % Sequential batch control 
Recipe_water_sp = [ 0  500  100  0  400  150  250   0    100 ];
for SQ = 1:1:length(Recipe_water_sp) % Fg Sequence number 
% User input aeration rates
if k<=Recipe_water(SQ)
    Fw =  Recipe_water_sp(SQ);
    break
else 
    Fw =  Recipe_water_sp(end);
end 
end
%% SBC - F_PAA
Recipe_PAA =    [25 200 1000 1500  1750 ]; % Sequential batch control 
Recipe_PAA_sp = [5  0    10  4  0  ];
for SQ = 1:1:length(Recipe_PAA_sp) % Fg Sequence number 
% User input aeration rates
if k<=Recipe_PAA(SQ)
    Fpaa =  Recipe_PAA_sp(SQ);
    break
else 
    Fpaa =  Recipe_PAA_sp(end);
end 
end
% Add PRBS to substrate flow rate  (Fpaa)
if Ctrl_flags.PRBS== 1;
if k > 500 && rem(k,100)==0
    rng('shuffle');
random_number = randi(3,1,1);
    noise_factor =1; 
if random_number == 1; 
    random_noise = 0;
elseif random_number == 2;
    random_noise = noise_factor;
else 
    random_noise = -noise_factor;
end 
X.PRBS_noise_addition(k) =random_noise;
end 
if k >475
     Fpaa = X.Fpaa.y(k-1);
end
if  k>500 && rem(k,100)==0
Fpaa = X.Fpaa.y(k-1) + X.PRBS_noise_addition(end);
end 
else 
    X.PRBS_noise_addition(k) =0;
    
end 



% SBC - NH3_shots 
Xd.NH3_shots.y(k)= 0;
end

%% Process faults
              % 0 - No Faults 
              % 1 - Aeration rate fault
              % 2 - Vessel back pressure  fault 
              % 3 - Substrate feed rate fault  
              % 4 - Base flowrate fault 
              % 5 - Coolant flowrate fault  
              % 6 - All of the above faults   
% Aeration fault
if Ctrl_flags.Faults ==1 ||  Ctrl_flags.Faults ==6
    if k>= 100 && k<=120
        Fg = 20; 
        u.Fault_ref = 1; % Fault occured
    end 
    if k>= 500 && k<= 550
        Fg = 20; 
        u.Fault_ref = 1; % Fault occured
    end 
    
end 
% Pressure fault
if Ctrl_flags.Faults ==2 ||  Ctrl_flags.Faults ==6
    if k>= 500 && k<=520
        pressure = 2; 
        u.Fault_ref = 1; % Fault occured
    end 
    if k>= 1000 && k<= 1200
        pressure = 2; 
        u.Fault_ref = 1; % Fault occured
    end 
end 
% Substrate feed fault
if Ctrl_flags.Faults ==3 ||  Ctrl_flags.Faults ==6
    if k>= 100 && k<=150
        Fs = 2; 
        u.Fault_ref = 1; % Fault occured
    end 
    if k>= 380 && k<= 460
        Fs = 20; 
        u.Fault_ref = 1; % Fault occured 
    end 
    if k>= 1000 && k<= 1070
        Fs = 20; 
        u.Fault_ref = 1; % Fault occured
    end 

% Process distrubance Fs shut down
% if k>= 500 && k<=1150
%         Fs = 0; 
%         u.Fault_ref = 1; % Fault occured
% end 
end 


% Base flowrate  fault
if Ctrl_flags.Faults ==4 ||  Ctrl_flags.Faults ==6
    if k>= 400 && k<=420
        Fb = 5; 
        u.Fault_ref = 1; % Fault occured
    end 
    if k>= 700 && k<= 800
        Fb = 10;
        u.Fault_ref = 1; % Fault occured
    end 
end 

% Coolant water  flowrate  fault
if Ctrl_flags.Faults ==5 ||  Ctrl_flags.Faults ==6
    if k>= 350 && k<=450
        Fc = 2; 
        u.Fault_ref = 1; % Fault occured
    end 
    if k>= 1200 && k<= 1350
        Fc = 10; 
        u.Fault_ref = 1; % Fault occured
    end 
end  

if Ctrl_flags.Raman_spec ==2
% Bulidng PID controller for PAA 
% builds the error history.  Samples 1 and 2 are calculated separately 
% because there is only instance of the error available
PAA_sp =  1200; 
if(k==1 || k==2)
    PAA_err  = PAA_sp-X.PAA.y(1);
    PAA_err1 = PAA_sp-X.PAA.y(1);
else
    PAA_err  = PAA_sp-X.PAA.y(k-1);
    PAA_err1 = PAA_sp-X.PAA.y(k-2);
end

% builds the temperature history of current and previous two samples.
if (k*h < 10)
    Fpaa = Fpaa;
    % after 80 hours use closed loop PI control
else
%     PAA_err  = PAA_sp-X.PAA_pred.y(k-3);
%     PAA_err1 = PAA_sp-X.PAA_pred.y(k-4);
if(k==1 || k==2)
    temp = X.PAA_pred.y(1);
    temp1 = X.PAA_pred.y(1);
    temp2 = X.PAA_pred.y(1);
elseif (k==3)
    temp = X.PAA_pred.y(2);
    temp1 = X.PAA_pred.y(1);
    temp2 = X.PAA_pred.y(1);
else
    temp = X.PAA_pred.y(k-2);
    temp1 = X.PAA_pred.y(k-3);
    temp2 = X.PAA_pred.y(k-4);
end
% PID control of PAA flowrate through Raman spectroscpy measurments
    if(k==1)    % 
        [Fpaa]= PIDSimple3(X.Fpaa.y(1), PAA_err, PAA_err1, temp, temp1, temp2, 0, 150, 0.1, 0.50, 0*0.002, h);   
    else
        [Fpaa] = PIDSimple3(X.Fpaa.y(k-1), PAA_err, PAA_err1, temp, temp1, temp2, 0, 150, 0.1, 0.50, 0*0.002, h);
    end
    
end
end


%% Controller vector
u.Fg = Fg;  % aeration rate [L/h]
u.RPM = 100;      % agitator RPM 
u.Fs = Fs;    % sugar feed rate [L/h]
u.Fa = Fa;          % acid flow rate [L/h]
u.Fb = Fb;          % base flow rate [L/h]
u.Fc = Fc;          % cooling water flow rate [L/h]
u.Fh = Fh;          % heating water flow [L/h]
u.d1 = ph_on_off;   % acid/base on off switch [0/1]
u.tfl = temp_on_off;% hot/cold water on off switch [0/1]
u.Fw = Fw;        % water for injection/dilution [L/h]
u.pressure  = pressure;  % air head pressure [bar] 
u.viscosity = viscosity;    % broth viscosity input [cPoise]
u.Fremoved = F_discharge;     % dumped broth flow [L/h]
u.Fpaa = Fpaa;               % PAA flow rate [L/h]
u.Foil = Foil;               % Foil flow rate [L/h]
u.NH3_shots = Xd.NH3_shots.y(k); % NH3 shots  added [kg]
% Defining Fault reference 
 if isfield(u, 'Fault_ref') ==1
     
 else 
     u.Fault_ref =0; 
 end 
     
%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037





