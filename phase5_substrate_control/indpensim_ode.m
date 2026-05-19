%% Industrial Penicillin Simulation indpensim_ode  
% This file is called by indpensim_run
% It contains the ordinary differential equations for the IndPenSim
% I/O structure : dy = indpensim_ode(t,y,FLAG,inp1);
%
%% Controller Inputs: 
% t: time span of integration 
% y: initial conditions for that time span 
% FLAG=[]; Matlab notation
% inp1: (1-by-26) vector of input parameters 
% inp1(1)  - Inhibition flag
% inp1(2)  - Fs Sugar feed rate (L/h)
% inp1(3)  - Fg Aeration rate  (m^(3)/h)
% inp1(4)  - RPM
% inp1(5)  - Fc Flow rate of cooling water (L/h)
% inp1(6)  - Fh Flow rate of heating water (L/h)
% inp1(7)  - Fb Flow rate of base (L/h)
% inp1(8)  - Fa Flow rate of acid (L/h)
% inp1(9)  - h_ode  ode solver step size (hours)
% inp1(10)  - Fw Flow rate of water for injection [L/h]
% inp1(11) - pressure Vessel head pressure (bar)
% inp1(12) - viscosity (recorded) [cP]
% inp1(13) - F_discharge Flow rate of discharge [L/h]
% inp1(14) - F_PAA Flow rate of Phenylacetic acid (L/hr)
% inp1(15) - F_oil Flow rate of soybean oil (L/hr)
% inp1(16) - F_NH3_shots Ammonia shots (kg)
% inp1(17) - Disturbance flag
% inp1(18) - dmuP Disturbance on penicillin growth rate [hr^{-1}] 
% inp1(19) - dmuX Disturbance on biomass growth rate [hr^{-1}]
% inp1(20) - distcs Disturbance on sugar inlet conc [g/L]
% inp1(21) - distcoil Disturbance on oil inlet conc [g/L]
% inp1(22) - distabc Disturbance on acid/base inlet conc [mol/L]
% inp1(23) - distPAA Disturbance on PAA inlet conc [mg/L]
% inp1(24) - distTcin Disturbance on coolant inlet temp [K]
% inp1(25) - distO_2in Disturbance on DO_2 inlet conc [mg/L]
% inp1(26) - Viscosity flag


%% ODE outputs:
% dy: matrix of time stamps and ODE solutions
%
%   Y(1) - X.S             substrate concentration [g/L]
%   Y(2) - X.DO2           dissolved oxygen concentration [mg/L]
%   Y(3) - X.O2            oxygen off gas []
%   Y(4) - X.P             penicillin concentration [g/L]
%   Y(5) - X.V             volume [L]
%   Y(6) - X.Wt            weight [kg]
%   Y(7) - X.pH            pH
%   Y(8) - X.T             temperature [K]
%   Y(9) - X.Q             generated heat [kcal]
%   Y(10) - X.viscosity    viscosity [cP]
%   Y(11) - A0+A2+A3+A4    Integral of biomass regions [g/L]
%   Y(12) - X.A0           Growing biomass region A0 [g/L]
%   Y(13) - X.A1           Non-growing biomass region A1 [g/L]
%   Y(14) - X.A3           Degenerated regions A3 [g/L]
%   Y(15) - X.A4           Autolysed regions A4 [g/L]
%   Y(16) -> Y(25)         mean vacoule number density functions (# cm^(-1) L)        
%   Y(26)  dn_m_dt         maximum vacuole volume department
%   Y(27)  dphi_0_dt       mean vacuole volume
%   Y(28) - CO_2           off-gas concentration [% off-gas ]
%   Y(29) - CO_2d          dissolved CO_2 [g/L]   
%   Y(30) - PAA            Phenlyacetic acid concentration [mg/L]
%   Y(31) - NH3            Nitrogen concentration [mg/L]
%   Y(32) - mup            Current penicillin growth rate [h^{-1]
%   Y(33) - mux            Current biomass growth rate [h^{-1]
function dy=indpensim_ode(t,y,asa,inp1,par)
mu_p                =       par(1);
mux_max             =       par(2);
ratio_mu_e_mu_b     =       par(3);
P_std_dev           =       par(4);
mean_P              =       par(5);
mu_v                =       par(6);
mu_a                =       par(7);
mu_diff             =       par(8);
beta_1              =       par(9);
K_b                 =       par(10);
K_diff              =       par(11);
K_diff_L            =       par(12);
K_e                 =       par(13);
K_v                 =       par(14);
delta_r             =       par(15);
k_v                 =       par(16);
D                   =       par(17);
rho_a0              =       par(18);
rho_d               =       par(19);
mu_h                =       par(20);
r_0                 =       par(21);
delta_0             =       par(22);

%% Process related parameters
Y_sX                =       par(23);
Y_sP                =       par(24);
m_s                 =       par(25);
c_oil               =       par(26);
c_s                 =       par(27);
Y_O2_X              =       par(28);
Y_O2_P              =       par(29);
m_O2_X              =       par(30);
alpha_kla           =       par(31);
a                   =       par(32);
b                   =       par(33);
c                   =       par(34);
d                   =       par(35);
Henrys_c            =       par(36);
n_imp               =       par(37);
r                   =       par(38);
r_imp               =       par(39);
Po                  =       par(40);
epsilon             =       par(41);
g                   =       par(42);
R                   =       par(43);
X_crit_DO2          =       par(44);
P_crit_DO2          =       par(45);
A_inhib             =       par(46);
Tf                  =       par(47);
Tw                  =       par(48);
Tcin                =       par(49);
Th                  =       par(50);
Tair                =       par(51);
C_ps                =       par(52);
C_pw                =       par(53);
dealta_H_evap       =       par(54);
U_jacket            =       par(55);
A_c                 =       par(56);
Eg                  =       par(57);
Ed                  =       par(58);
k_g                 =       par(59);
k_d                 =       par(60);
Y_QX                =       par(61);
abc                 =       par(62);
gamma1              =       par(63);
gamma2              =       par(64);
m_ph                =       par(65);
K1                  =       par(66);
K2                  =       par(67);
N_conc_oil          =       par(68);
N_conc_paa          =       par(69);
N_conc_shot         =       par(70);
Y_NX                =       par(71);
Y_NP                =       par(72);
m_N                 =       par(73);
X_crit_N            =       par(74);
PAA_c               =       par(75);
Y_PAA_P             =       par(76);
Y_PAA_X             =       par(77);
m_PAA               =       par(78);
X_crit_PAA          =       par(79);
P_crit_PAA          =       par(80);
B_1                 =       par(81);
B_2                 =       par(82);
B_3                 =       par(83);
B_4                 =       par(84);
B_5                 =       par(85);
delta_c_0           =       par(86);
k3                 =       par(87);
k1                  =       par(88);
k2                  =       par(89);
t1                  =       par(90);
t2                  =       par(91);
q_co2               =       par(92);
X_crit_CO2          =       par(93);
alpha_evp           =       par(94);
beta_T              =       par(95);
pho_g               =       par(96);
pho_oil             =       par(97);
pho_w               =       par(98);
pho_paa             =       par(99);
O_2_in              =       par(100);
N2_in               =       par(101);
C_CO2_in            =       par(102);
Tv                  =       par(103);
T0                  =       par(104);
alpha_1             =       par(105);

%% Process inputs
inhib_flag=inp1(1); % Inhibition flag
Fs=inp1(2); % Substrate feed rate [L/h]
Fg=(inp1(3)/60); % Aeration rate [m^3 s^{-1}]
RPM=inp1(4); % RPM 
%Tcin=inp1(5); % Temperature of cooling water in [K]
Fc=inp1(5); % cooling water jacket flow rate [L/h]
Fh = inp1(6); % heating water jacket flow rate [L/h]
Fb=inp1(7); % Base flow rate [L/h]
Fa=inp1(8); % Acid flow rate [L/h]
%abc=inp1(10); % Acid/base concentration [moles]
step1 = inp1(9); % ode step size [h]
Fw = inp1(10); % water for injection [L/h]
Fw(Fw<0) = 0; % correcting for negative flow values
pressure = inp1(11); % vessel head pressure [bar]
% Viscosity flag 
if inp1(26) == 0;% Ctrl_flags.Vis= 0 - Uses simulated viscosity 
viscosity = y(10);    % viscosity [cP] 
else %Ctrl_flags.Vis= 1 - Uses viscosity recorded in batch records as input
viscosity = inp1(12); % viscosity [cP]
end 
F_discharge = inp1(13); % discharge flow rate [L/h]
Fpaa = inp1(14); % Phenylacetic acid flow   [L/h]
Foil = inp1(15); % Oil flow rate  [L/h]
NH3_shots = inp1(16); % Ammonia sulphate shots added [kg]
dist_flag = inp1(17); % Disturbance on/off flag [-]
distMuP =  inp1(18); % Disturbance on penicillin growth rate [hr^{-1}]
distMuX =  inp1(19); % Disturbance on biomass growth rate [hr^{-1}] 
distsc =  inp1(20); % Disturbance on sugar concentration [g L^{-1}] 
distcoil =  inp1(21); % Disturbance on oil concentration [g L^{-1}] 
distabc = inp1(22); % Disturbance on acid/base inlet concentration  [g L^{-1}]
distPAA = inp1(23); % Disturbance on PAA inlet concentration  [g L^{-1}]
distTcin = inp1(24); % Disturbance on inlet coolant temperature [K]
distO_2_in = inp1(25); % Disturbance on inlet oxygen concentration [%]
pho_b = (1100 + y(4) + y(12) + y(13) + y(14) + y(15)) ; % Broth viscosity
%Ctrl_flags.Dis = 0; (inp1(20)) % 0 - No disturbances 1- disturbances
if dist_flag == 1
mu_p    = mu_p + distMuP; % Penicillin specific growth rate + disturbance [hr^{-1}]
mux_max =  mux_max+ distMuX; %Biomass specific growth rate + disturbance hr^{-1}
c_s=c_s + distsc ; % Substrate conc [g/L]
c_oil = c_oil + distcoil; 
abc = abc +distabc;
PAA_c = PAA_c+ distPAA;
Tcin = Tcin + distTcin;
O_2_in = O_2_in+ distO_2_in;
end


%% Process parameters
% Adding in age-dependant term
A_t1 = ((y(11))/(y(12)+y(13)+y(14)+y(15))); %*factor_1*4;
% Variables
s = y(1); % substrate  g L^(-1)
a_1 = y(13); % Biomass (Extension) region -  A_1 [g L^{-l}]
a_0 = y(12); % Biomass (Branching) region - A_0 [g L^{-l}]
a_3 =y(14);  % Biomass 
total_X = y(12) + y(13) +y(14)+y(15); % Total Biomass 

% Calculating liquid height in vessel
 h_b = (y(5)/1000)/(pi()*(r^2));
 h_b = h_b*(1-epsilon); % (ungassed height)
% Calculating log mean pressure of vessel [bar]
pressure_bottom  =  1+ pressure + ((pho_b*h_b)*9.81*10^(-5)); % [bar]
pressure_top = 1+ pressure; % [bar]
log_mean_pressure = (pressure_bottom - pressure_top)./(log(pressure_bottom./pressure_top));
total_pressure = log_mean_pressure;
% Ensuring minimum value for viscosity
viscosity(viscosity<4) = 1;
DOstar_tp =  (((total_pressure)*O_2_in)/Henrys_c);%in mg/l (Henry's constant has the units of (bar L mg-1)

%% Inhibition flags 
if inhib_flag == 0 
  pH_inhib = 1; 
  NH3_inhib = 1; 
  T_inhib = 1;
  mu_h = 0.003;  
  DO_2_inhib_X =1; 
  DO_2_inhib_P =1; 
  CO2_inhib =1; 
  PAA_inhib_X= 1; 
  PAA_inhib_P =1; 
end 

if inhib_flag == 1 % Inhibition DO2, T, pH 
% adding in pH effect
pH_inhib =  (1/(1+ (y(7)/K1)+ (K2/y(7))));
% Ammonia inhibition on Biomass 
NH3_inhib = 1;
% Temperature inhibition    
T_inhib =  (k_g*exp(-(Eg/(R*y(8)))) - k_d*exp(-(Ed/(R*y(8)))))*0+1;
% Carbon Dioxide inhibition
CO2_inhib =  1;

 %Dissolved oxygen inhibition
DO_2_inhib_X =0.5*(1-tanh(A_inhib*(X_crit_DO2*(((total_pressure)*O_2_in)/Henrys_c)-y(2)))); 
DO_2_inhib_P =0.5*(1-tanh(A_inhib*(P_crit_DO2*(((total_pressure)*O_2_in)/Henrys_c)-y(2)))); 
% Phenylacetic acid (PAA) inhibition terms
PAA_inhib_X = 1;
PAA_inhib_P = 1;

 % Temperature and pH effect on hydrolysis rate
pH = -log10(y(7));
k4 = exp((B_1 +B_2*pH+B_3*y(8)+ B_4*(pH^2))+B_5*(y(8)^2));
mu_h = k4;
end


if inhib_flag == 2 %Inhibition of DO2,T,pH,CO_2_L,PAA and N  
% adding in pH effect
pH_inhib =  1/(1+ (y(7)/K1)+ (K2/y(7)));
% Ammonia inhibition on Biomass 
NH3_inhib = 0.5*(1-tanh(A_inhib*(X_crit_N-y(31))));
% Temperature inhibition    
T_inhib =  k_g*exp(-(Eg/(R*y(8)))) - k_d*exp(-(Ed/(R*y(8)))) ;
% Carbon Dioxide inhibition
CO2_inhib =   0.5*(1+tanh(A_inhib*(X_crit_CO2-y(29)*1000)));
 %Dissolved oxyegn inhibition
DO_2_inhib_X = 0.5*(1-tanh(A_inhib*(X_crit_DO2*(((total_pressure)*O_2_in)/Henrys_c)-y(2)))); 
DO_2_inhib_P = 0.5*(1-tanh(A_inhib*(P_crit_DO2*(((total_pressure)*O_2_in)/Henrys_c)-y(2)))); 
% Phenylacetic acid (PAA) inhibition terms
PAA_inhib_X = 0.5*(1+(tanh((X_crit_PAA - y(30)))));
PAA_inhib_P = 0.5*(1+(tanh((-P_crit_PAA + y(30)))));
% Temperature and pH effect on hydrolysis rate
pH = -log10(y(7));
k4 = exp((B_1 +B_2*pH+B_3*y(8)+ B_4*(pH^2))+B_5*(y(8)^2));
mu_h = k4;
end





%% Main rate equations for kinetic expressions
% Penicillin inhibition curve
P_inhib = 2.5*P_std_dev*((P_std_dev*sqrt(2*pi))^-1*exp(-.5*((s-mean_P)/P_std_dev).^2));
% Specific growth rates of biomass regions with inhibition effect
mu_a0 = ratio_mu_e_mu_b*mux_max*pH_inhib*NH3_inhib*T_inhib*DO_2_inhib_X*CO2_inhib*PAA_inhib_X; % Rate constant for Branching A0
mu_e =                  mux_max*pH_inhib*NH3_inhib*T_inhib*DO_2_inhib_X*CO2_inhib*PAA_inhib_X; % Rate constant for extension A1

 K_diff = par(11)-(A_t1*beta_1);
if K_diff < K_diff_L
    K_diff = K_diff_L;
end
% Growing A_0 region
r_b0 =  mu_a0*a_1*s/(K_b+s); % rate equation 
r_sb0  = Y_sX*r_b0; % rate equation of substrate consumption  

% Non-growing regions A_1 region 
r_e1 = (mu_e*a_0*s)/(K_e+s); %rate equation 
r_se1 = Y_sX*r_e1; % rate equation of substrate consumption

% Differentiation (A_0 -> A_1)
r_d1 = mu_diff*a_0/(K_diff+s);
r_m0 =  m_s*a_0/(K_diff+s);
% rho_a0 = 0.35; % Vacuole density g cm^-3


n=17;
phi(1) = y(27);
for k= 2:1:10
r_mean(k)= (1.5e-4) + (k-2)*delta_r;
phi(k) = ((4*pi*r_mean(k)^3)/3)*y(n)*delta_r;
n=n+1;
end
% Total vacuole volume
v_2 = sum((phi));
rho_a1 = (a_1/((a_1/rho_a0)+ v_2)); % Density of non-growing regions (A_{1}) g mL^-1
v_a1 = a_1/(2*rho_a1) -v_2; % Total volume of non-growing regions
% Penicillin produced from the non-growing regions  A_1 regions  
r_p = mu_p*rho_a0*v_a1*P_inhib*DO_2_inhib_P*PAA_inhib_P - mu_h*y(4);


% ----- Vacuole formation-------
r_m1 = (m_s*rho_a0*v_a1*s)/(K_v +s);

% ------ Vacuole degeneration -------------------
r_d4 = mu_a*a_3; % Biomass autolysis rate

% ------ Vacuole Volume -------------------
% n_0 - mean vacoule number density for vacuoles sized ranging from delta_0
% -> r_0
dn0_dt = ((mu_v*v_a1)/(K_v +s))*((6/pi)*((r_0 +delta_0)^-3)) -k_v*y(16);  % y(16)

n = 17;
% n_j - mean vacoule number density for vacuoles sized ranging from r_{j}
% -> r_{j+1} where j = 1,2...9.
dn1_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2; %y(17)
n = n +1;
dn2_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2; %y(18)
n = n +1;
dn3_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2;%y(19)
n = n +1;
dn4_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2;%y(20)
n = n +1;
dn5_dt =-k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2;%y(21)
n = n +1;
dn6_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2; %y(22)
n = n +1;
dn7_dt =-k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2;%y(23)
n = n +1;
dn8_dt =-k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2; %y(24)
n = n +1;
dn9_dt = -k_v*((y(n+1)-y(n-1))/(2*(delta_r))) + D*(y(n+1)-2*y(n)+y(n-1))/(delta_r)^2; %y(25)
n_k = dn9_dt;
k = 10;
% Mean vacoule density for  department k all vacuoles above k in size are
% are assumed constant size
r_k=  r_0 + (k-2)*delta_r;
k = 12;
r_m= ( r_0+ (k-2)*delta_r);
% Calculating maximum vacuole volume department
dn_m_dt = k_v*n_k/(r_m-r_k) -mu_a*y(26); % y(26)
n_k = y(25);
% mean vacuole
dphi_0_dt = ((mu_v*v_a1)/(K_v +s)) -k_v*y(16)*(pi*(r_0 +delta_0)^3)/6;

%% Volume and Weight expressions
% Working out the Volume (V) and Vessel Weight (Wt)  
F_evp = y(5)*alpha_evp*(exp(2.5*(y(8) - T0)/(Tv-T0))-1); % Evaporation rate [L/hr]

pho_feed = (c_s/1000*pho_g +(1-c_s/1000)*pho_w); % [g/L^3]
% Dilution term [L]
dilution = Fs+ Fb+Fa+Fw-F_evp+ Fpaa;
% Change in Volume [L]
dV1 = Fs+ Fb+Fa+Fw+F_discharge/(pho_b/1000) -F_evp+ Fpaa; 
% Change in Weight [kg]
dWt  = Fs*pho_feed/1000+ pho_oil/1000*Foil+ Fb+Fa+Fw+F_discharge -F_evp + Fpaa*pho_paa/1000; 

%% ODE's for Biomass regions
% Growing regions  a_{0} [g/l]
da_0_dt = r_b0 - r_d1- y(12)*dilution/y(5);
% Non growing regions a_{1} [g/l]
da_1_dt  = r_e1 - r_b0+ r_d1 - (pi*((r_k +r_m)^3)/6)*rho_d*k_v*n_k - y(13)*dilution/y(5);
% Degenerated regions a_{3} [g/l]
da_3_dt = (pi*((r_k +r_m)^3)/6)*rho_d*k_v*n_k - r_d4 - y(14)*dilution/y(5);
% Autolysed regions a_{4} [g/l]
da_4_dt =  r_d4- y(15)*dilution/y(5);
% Penicillin production 
dP_dt = r_p -y(4)*dilution/y(5);

% Active Biomass rate [g/L hr ]
X_1 = da_0_dt + da_1_dt + da_3_dt+da_4_dt; 
% % Total biomass [g/L]
X_t = y(12)+y(13)+y(14)+y(15);


Qrxn_X = X_1*Y_QX*y(5)*Y_O2_X/1000;
Qrxn_P = dP_dt*Y_QX*y(5)*Y_O2_P/1000;
    
 Qrxn_t =  Qrxn_X+Qrxn_P;
% 
if Qrxn_t<0
    Qrxn_t =0;
end 
Qrxn_t = Qrxn_t;


N = RPM/60; % rotations per second
D_imp= 2*r_imp;
unaerated_power = (n_imp*Po*pho_b*(N^3)*(D_imp^5)); % 
P_g = 0.706*(((unaerated_power^2)*(N)*D_imp^3)/(Fg^0.56))^0.45;
P_n = P_g/unaerated_power;
variable_power = (n_imp*Po*pho_b*(N^3)*(D_imp^5)*P_n)/1000; % divide by 1000 for kW 

%------------------------------------------------------------------------------------------------------------
%                           Process parameters
%------------------------------------------------------------------------------------------------------------
%%---------ODEs of main variables 1-32

% Defining the column vector
dy=zeros(31,1);

%% Substrate utilization [g/L] 
dy(1) =-r_se1-r_sb0-r_m0-r_m1 -((Y_sP*mu_p*rho_a0*v_a1*P_inhib*DO_2_inhib_P*PAA_inhib_P))+Fs*c_s/y(5) + Foil*c_oil/y(5) -y(1)*dilution/y(5);

%% Dissolved oxygen [mg/L]
V_s = (Fg)/(pi*(r^2)); % superficial gas velocity [m/s]
T = y(8); % temperature [K]  
V = y(5); % Volume [V]
V_m = y(5)/1000; % Volume [m^{3}]
 P_air = ((V_s*R*T*V_m/(22.4*h_b))*log(1+pho_b*9.81*h_b/(pressure_top*10^5))); % aeration power rate [kW]
P_t1 = (variable_power+P_air); 
viscosity(viscosity<=4) = 1;
vis_scaled =(viscosity/100);% Scaled viscosity [cP]
 oil_f = (Foil/V);
 kla = (alpha_kla*(((V_s.^(a)).*((P_t1./V_m).^ b).*(vis_scaled).^(c))).*(1-oil_f.^(d)));
 OUR = (-X_1)*Y_O2_X-m_O2_X*X_t-dP_dt*Y_O2_P;
  OTR = kla*(DOstar_tp-y(2));
 dy(2) = OUR+OTR -(y(2)*dilution/y(5));

%% O_2 off-gas [%]  
Vg = epsilon*V_m; 
Qfg_in = 60*Fg*1000*32/22.4;   
Qfg_out = 60*Fg*(N2_in/(1-y(3) - y(28)/100))*1000*32/22.4;
dy(3) = (Qfg_in*O_2_in - Qfg_out*y(3) -0.001*OTR*V_m*60)/(Vg*28.97*1000/22.4); % Oxygen transfer rate g/L

%% Penicillin production rate [g/L]
dy(4) = r_p -y(4)*dilution/y(5);
%% Volume change [L]
dy(5)=dV1;
%% Weight change [Kg]
dy(6) = dWt;
%% pH 
% pH disturbances
pH_dis = Fs + Foil +Fb +Fa +F_discharge+Fw; %Removing water addition as doesn't effet the pH
% conversion between [H+] balance and [OH-]
if -log10(y(7)) < 7 % acidic = [H+] balance 
    cb = -abc; 
    ca = +abc; 
     y(7) = y(7);
    pH_balance = 0; 
else  % basic = [OH-] balance 
    cb = +abc; 
%     cb = 0.2;  
    ca = -abc;
    y(7) = (1e-14/y(7)-y(7)); 
    pH_balance = 1;      
end 
% Calculation of ion addition
  B=(y(7)*y(5)+ca*Fa*step1+cb*Fb*step1)/(y(5)+Fb*step1+Fa*step1);
  B = -B;
if pH_balance == 1; %basic
      %  dy(7)=-ck*(r_b0 +r_e1+r_d4+ r_d1+m_ph*total_X) - ck1*r_p - ck2*(pH_dis)+[(-B-sqrt(B^2+4e-14))/2-y(7)]/step1;
    dy(7)=-gamma1*(r_b0 +r_e1+r_d4+ r_d1+m_ph*total_X) - gamma1*r_p - gamma2*(pH_dis)+((-B-sqrt(B^2+4e-14))/2-y(7));
   end 
if pH_balance == 0; % acidic
   % dy(7)=+ck*(r_b0 +r_e1+r_d4+ r_d1+ m_ph*total_X) + ck1*r_p + ck2*(pH_dis)+[(-B+sqrt(B^2+4e-14))/2-y(7)]/step1;
     dy(7)=+gamma1*(r_b0 +r_e1+r_d4+ r_d1+m_ph*total_X) + gamma1*r_p + gamma2*(pH_dis)+((-B+sqrt(B^2+4e-14))/2-y(7));
      % dy(7)=(1e-14/y(7)-y(7)); 
end 

%% Temperature  [K]
Ws = P_t1;
Qcon = U_jacket*A_c*(y(8)-Tair); 
dQ_dt = Fs*pho_feed*C_ps*(Tf-y(8))/1000 +Fw*pho_w*C_pw*(Tw-y(8))/1000 ...
    - F_evp*pho_b*C_pw/1000 - dealta_H_evap*F_evp*pho_w/1000 ...
    +Qrxn_t+Ws-(alpha_1/1000)*Fc^(beta_T+1)...
     *((y(8)-Tcin)/(Fc/1000+(alpha_1*(Fc/1000)^beta_T)/2*pho_b*C_ps)) -(alpha_1/1000)*Fh^(beta_T+1)...
     *((y(8)-Th)/(Fh/1000+(alpha_1*(Fh/1000)^beta_T)/2*pho_b*C_ps))   -Qcon;

dy(8) = dQ_dt/((y(5)/1000)*C_pw*pho_b) ; % Temperature 

% Heat generation [kJ]   
dy(9)=dQ_dt;
%% Viscosity 

dy(10) =  (((3*(a_0^(1/3))*(1/(1+ exp(-k1*(t-t1))))*(1/(1+ exp(-k2*(t-t2)))) - k3*Fw)));

%% Total X [g/L]
dy(11) = (y(12)+y(13)+y(14)+y(15)); 
%%--------------------------------------------------------------------------
%         Adding in the ODE's for hyphae
%--------------------------------------------------------------------------
% da0/dt - y(12)  - Growing regions 
dy(12) = da_0_dt;
% da1/dt - y(13)  - Non -growing regions
dy(13)  = da_1_dt;
% da3/dt - y(14)  - degenerated regions
dy(14)  = da_3_dt;
% da4/dt autolysed biomass
dy(15) =  da_4_dt;
% Vacuole regions
dy(16) = dn0_dt;
dy(17) = dn1_dt;
dy(18) = dn2_dt;
dy(19) = dn3_dt;
dy(20) = dn4_dt;
dy(21) = dn5_dt;
dy(22) = dn6_dt;
dy(23) = dn7_dt;
dy(24) = dn8_dt;
dy(25) = dn9_dt;
dy(26) = dn_m_dt;
dy(27) = dphi_0_dt;
%% CO_2 [%]
total_X_CO2 = y(12) + y(13);
CER = total_X_CO2*q_co2*V;
dy(28) =((((60*Fg*44*1000)/22.4)*C_CO2_in + CER - ((60*Fg*44*1000)/22.4)*y(28)))/(Vg*28.97*1000/22.4);
% %% dissolved CO_2 [mg/L]
% Henrys_c_co2 = (exp(11.25 - 395.9   /(y(8) - 175.9))/44)/100; % convert to g CO_2/ L bar
% C_star_CO2 = total_pressure*Henrys_c_co2*y(28)/100;
% dy(29) = kla*delta_c_0*(C_star_CO2-y(29)) - y(29)*dilution/y(5);
%% dissolved CO_2 [mg/L]
Henrys_c_co2 = (exp(11.25 - 395.9/(y(8) - 175.9)))/(44*100); % convert to g CO_2/ L bar
C_star_CO2 = (total_pressure*y(28))/Henrys_c_co2;
dy(29) = kla*(delta_c_0)*(C_star_CO2-y(29)) - y(29)*dilution/y(5);


%% PAA [mg/L]
dy(30) =Fpaa*PAA_c/V   -(Y_PAA_P*dP_dt) -Y_PAA_X*X_1  -m_PAA*y(4)   -y(30)*dilution/y(5); %PAA
%% N [mg/L]
X_C_nitrogen = (-r_b0-r_e1-r_d1 -r_d4)*Y_NX;
P_C_nitrogen = -dP_dt*Y_NP;
dy(31) = (NH3_shots*N_conc_shot)/y(5) +X_C_nitrogen + P_C_nitrogen - m_N*total_X +(1*N_conc_paa*Fpaa/y(5)) +N_conc_oil*Foil/y(5) - y(31)*dilution/y(5);
dy(32) = mu_p;
dy(33) = mu_e;
end
%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037
