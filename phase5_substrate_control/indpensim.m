%% Penicillin G production simulator indpensim.m
%  Usage:
%  X = indpensim(f_input, Xd, x0, h, T)
%  Parameters:
%  u - @f_input(X, Xd, k, h, T) function that at each sampling instance (defined
%      by h and T) returns the values of the manipulated variables.
%      X is the current batch data. Xd is the industrial data(if imported). See next
%      section for the structure of u.
%
% 
%  x0 - structure containing the initial states. See next section for
%       details.
%  h -  sampling time (hours)
%  T -  experiment length, must be multiple of sampling time
%  X -  time history of variables used in the simulation. See the next
%       section for details.
%
%  Input, state and output variables description
%   x0       structure with initial conditions
%   x0.S     Initial substrate concentration [g/L]
%   x0.DO2   Initial dissolved oxygen concentration [mg/L]
%   x0.O2    Initial oxygen off gas [%]
%   x0.P     Initial penicillin concentration [g/L]
%   x0.V     Initial volume [L]
%   x0.WT    Initial weight [kg]
%   x0.pH    Initial pH
%   x0.T     Initial temperature [K]
%   x0.a0    Initial A0 biomass concentration [g/L]
%   x0.a1    Initial A1 biomass concentration [g/L]
%   x0.a3    Initial A3 biomass concentration [g/L]
%   x0.a4    Initial A4 biomass concentration [g/L]

%   h - simulation sampling period
%   T - simulation length
%   X - batch structure with the simulation results

%   IndPenSim Manipulated Variables:
%   X.Fg            aeration rate [m^{3}/min]
%   X.Pw            aeration power [kW]
%   X.Fs            substrate feed rate [L/h]
%   X.Foil          soybean oil feed rate [L/h]
%   X.Tf            substrate feed temperature [K]
%   X.Fa            acid flow rate [L/h]
%   X.Fb            base flow rate [L/h]
%   X.Fc            heating/cooling water flow rate [L/h]
%   X.RPM           Agitator RPM
%   IndPenSim Manipulated Variables:
%   X.Fw        water for injection/dilution [L/h]
%   X.pressure      vessel pressure [bar] 
%   X.viscosity     broth viscosity input [cPoise]
%   X.Fremoved      dumped broth flow [L/h]

%   IndPenSim States: [Y(1) -> Y(34)]
%   Y(1) - X.S             substrate concentration [g/L]
%   Y(2) - X.DO2           dissolved oxygen concentration [mg/L]
%   Y(3) - X.O_2           O_2 off-gas concentration [% off-gas ]
%   Y(4) - X.P             penicillin concentration [g/L]
%   Y(5) - X.V             volume [L]
%   Y(6) - X.Wt            weight [kg]
%   Y(7) - X.pH            pH
%   Y(8) - X.T             temperature [K]
%   Y(9) - X.Q             generated heat [kJ]
%   Y(10) - X.viscosity    viscosity [cP]
%   Y(11) - A0+A2+A3+A4    Integral of biomass regions [g/L h]
%   Y(12) - X.A0           Growing biomass region A0 [g/L]
%   Y(13) - X.A1           Non-growing biomass region A1 [g/L ]
%   Y(14) - X.A3           Degenerated regions A3 [g/L]
%   Y(15) - X.A4           Autolysed regions A4 [g/L]
%   Y(16) -> Y(25)         mean vacoule number density function (# cm^(-1) L)        
%   Y(26)  dn_m_dt         maximum vacuole volume department
%   Y(27)  dphi_0_dt       mean vacuole volume
%   Y(28) - CO_2           off-gas concentration [% off-gas ]
%   Y(29) - CO_2d          dissolved CO_2 [mg/L]   
%   Y(30) - PAA            Phenlyacetic acid concentration [mg/L]
%   Y(31) - NH3            Nitrogen concentration [mg/L]
%   IndPenSim states:
% u - structure with the values of the manipulated variables
%   u.Fg            aeration rate [L/s]
%   u.Pw            aeration power [kW]
%   u.Fs            sugar feed rate [L/h]
%   u.Tf            substrate feed temperature [K]
%   u.Fa            acid flow rate [L/h]
%   u.Fb            base flow rate [L/h]
%   u.Fc            heating/cooling water flow rate [L/h]
%   u.d1            acid/base on off switch [0/1]
%   u.Fw            water for injection/dilution [L/h]
%   u.pressure      vessel head pressure [bar] 
%   u.viscosity     broth viscosity input [cPoise]
%   u.Fremoved      Discharge rate [L/h]


%% indpensim code
function [X] = indpensim(f_input, Xd, x0, h, T, solv, p, Ctrl_flags)

%% simulation timing init
N = T/h; %experiment length in samples
h_ode = h/20; % ode solver step size (hours)
t = 0:h:T; % time vector

% creates batch structure 
X = createBatch(h,T);

% User control inputs; 
% converts from pH to H+ conc.
x0.pH=10^(-x0.pH);

%% Main ode loop  
for k=1:1:N
%    waitbar(k / N);
    % fills the batch with just the initial conditions so the control system
    % can provide the first input. These will be overwritten after
    % the ODEs are integrated.
    if(k==1)
        X.S.y(1)   = x0.S;
        X.DO2.y(1) = x0.DO2;
        X.X.y(1)   = x0.X;
        X.P.y(1)   = x0.P;
        X.V.y(1)   = x0.V;
        X.CO2outgas.y(1) = x0.CO2outgas;
       % X.CO2.y(1) = x0.CO2;
        X.pH.y(1)  = x0.pH;
        X.T.y(1)   = x0.T;
    end
    
    % gets MVs
    [u X] = f_input(X, Xd, k, h, T,Ctrl_flags);
    
    % builds initial conditions and control vectors specific to
    % indpensim_ode using ode45
          if(k==1)
             x00 = [x0.S x0.DO2 x0.O2 x0.P x0.V x0.Wt x0.pH x0.T 0 4 x0.Culture_age x0.a0 x0.a1 x0.a3 x0.a4 zeros(1,12) x0.CO2outgas 0 x0.PAA x0.NH3 0 0 ];
         else
            x00 = [ X.S.y(k-1)...
                X.DO2.y(k-1)...
                X.O2.y(k-1) ...
                X.P.y(k-1)...
                X.V.y(k-1)...
                X.Wt.y(k-1)... 
                X.pH.y(k-1)...
                X.T.y(k-1)...
                X.Q.y(k-1)... 
                X.Viscosity.y(k-1)... 
                X.Culture_age.y(k-1)...
                X.a0.y(k-1)... 
                X.a1.y(k-1)...
                X.a3.y(k-1)...
                X.a4.y(k-1)...
                X.n0.y(k-1)...
                X.n1.y(k-1)...
                X.n2.y(k-1)...
                X.n3.y(k-1)...
                X.n4.y(k-1)...
                X.n5.y(k-1)...
                X.n6.y(k-1)...
                X.n7.y(k-1)...
                X.n8.y(k-1)...
                X.n9.y(k-1)...
                X.nm.y(k-1)...
                X.phi0.y(k-1)...
                X.CO2outgas.y(k-1)...
                X.CO2_d.y(k-1) ...
                X.PAA.y(k-1) ...
                X.NH3.y(k-1) ...
                0 ...
                0];
        end
        
        %% Process disturbances
        distMuP = Xd.distMuP.y(k);
        distMuX = Xd.distMuX.y(k);
        distcs  = Xd.distcs.y(k);
        distcoil = Xd.distcoil.y(k);
        distabc = Xd.distabc.y(k);
        distPAA = Xd.distPAA.y(k);
        distTcin = Xd.distTcin.y(k);
        distO_2in = Xd.distO_2in.y(k);
        
        u00 = [Ctrl_flags.Inhib u.Fs u.Fg u.RPM  u.Fc u.Fh u.Fb u.Fa  h_ode u.Fw u.pressure u.viscosity   u.Fremoved u.Fpaa u.Foil  u.NH3_shots   Ctrl_flags.Dis distMuP distMuX distcs distcoil distabc distPAA distTcin distO_2in Ctrl_flags.Vis];
% To account for inability of growth rates of biomass and penicillin to 
%return to normal after continuous periods of suboptimal pH and temperature conditions
%If the Temperature or pH results is off set-point for k> 100 mu_p(max) is reduced to current value 
 if Ctrl_flags.Inhib ==1 || Ctrl_flags.Inhib ==2  
  % If the Temperature or pH results is off set-point for k> 40 mu_x(max) is reduced to current value 
% if k > 40
%         a1 =  diff((X.mu_X_calc.y((k-40:k-1),1)));
%         a2 = a1(2:end)>0;
%           if sum(a2) <= 1    
%          p(2) = X.mu_X_calc.y(k-1).*5;
%           end 
% end    

 if k > 65
      a1 =  diff((X.mu_X_calc.y((k-65:k-1),1)));
      a2 = a1<0;
      if sum(a2) >= 63
          p(2) = X.mu_X_calc.y(k-1).*5; % Changing mu_X to current minimum value
      end
  end

    
end
    
        
        %% Solver selection and calling indpensim_ode
    if(solv==1)
     [t_sol,y_sol] = ode45('indpensim_ode', t(k):h_ode:t(k+1), x00, [], u00,p);
    end 
     if(solv==2)
          [t_sol,y_sol] = ode15s('indpensim_ode', t(k):h_ode:t(k+1), x00,[], u00, p);
             warning off; % disables warning message for integration tolerance
     end
     if(solv==3)
         [t_sol,y_sol] = ode23t('indpensim_ode', t(k):h_ode:t(k+1), x00, [], u00,p);
     end
      % Defining minimum value for all variables for numerical stability
    for n =1:1:31
    if y_sol(end,n)<=0
        y_sol(end,n)=0.001;
    end 
    end 

    %% Saving all  manipulated variables
    X.Fg.t(k) = t_sol(end);
    X.Fg.y(k) = u.Fg;
    X.RPM.t(k) = t_sol(end);
    X.RPM.y(k) = u.RPM;
    X.Fpaa.t(k) = t_sol(end);
    X.Fpaa.y(k) = u.Fpaa;
    X.Fs.t(k)  = t_sol(end);
    X.Fs.y(k)  = u.Fs;
    X.Fa.t(k) = t_sol(end);
    X.Fa.y(k) = u.Fa;
    X.Fb.t(k) = t_sol(end);
    X.Fb.y(k) = u.Fb;
    X.Fc.t(k) = t_sol(end);
    X.Fc.y(k) = u.Fc;
    X.Foil.t(k) = t_sol(end);
    X.Foil.y(k) = u.Foil;
    X.Fh.t(k) = t_sol(end);
    X.Fh.y(k) = u.Fh;
    X.Fw.t(k) = t_sol(end);
    X.Fw.y(k) = u.Fw;
    X.pressure.t(k) = t_sol(end);
    X.pressure.y(k) = u.pressure;
    X.Fremoved.t(k) = t_sol(end);
    X.Fremoved.y(k) = u.Fremoved;
    
    %% Saving all the  IndPenSim states
    X.S.y(k)   = y_sol(end,1);
    X.S.t(k)   = t_sol(end);
    X.DO2.y(k) = y_sol(end,2);
  %  Required for numerical stability
    if    X.DO2.y(k)<2
         X.DO2.y(k) = 1; 
    else
          X.DO2.y(k) =   X.DO2.y(k);
    end 
    X.DO2.t(k) = t_sol(end);
    X.O2.y(k)   = y_sol(end,3);
    X.O2.t(k)   = t_sol(end);
    X.P.y(k)   = y_sol(end,4);
    X.P.t(k)   = t_sol(end);
    X.V.y(k)   = y_sol(end,5);
    X.V.t(k)   = t_sol(end);
    X.Wt.y(k) = y_sol(end,6);
    X.Wt.t(k) = t_sol(end);
    X.pH.y(k)  = y_sol(end,7);
    X.pH.t(k)  = t_sol(end);
    X.T.y(k)   = y_sol(end,8);
    X.T.t(k)   = t_sol(end);
    X.Q.y(k)   = y_sol(end,9);
    X.Q.t(k)   = t_sol(end);
    X.Viscosity.y(k)   = y_sol(end,10);
    X.Viscosity.t(k)   = t_sol(end);
    X.Culture_age.y(k)   = y_sol(end,11);
    X.Culture_age.t(k)   = t_sol(end);
    X.a0.y(k)   = y_sol(end,12);
    X.a0.t(k)   = t_sol(end);
    X.a1.y(k)   = y_sol(end,13);
    X.a1.t(k)   = t_sol(end);
    X.a3.y(k)   = y_sol(end,14);
    X.a3.t(k)   = t_sol(end);
    X.a4.y(k)   = y_sol(end,15);
    X.a4.t(k)   = t_sol(end);
    X.n0.y(k)   = y_sol(end,16);
    X.n0.t(k)   = t_sol(end);
    X.n1.y(k)   = y_sol(end,17);
    X.n1.t(k)   = t_sol(end);
    X.n2.y(k)   = y_sol(end,18);
    X.n2.t(k)   = t_sol(end);
    X.n3.y(k)   = y_sol(end,19);
    X.n3.t(k)   = t_sol(end);
    X.n4.y(k)   = y_sol(end,20);
    X.n4.t(k)   = t_sol(end);
    X.n5.y(k)   = y_sol(end,21);
    X.n5.t(k)   = t_sol(end);
    X.n6.y(k)   = y_sol(end,22);
    X.n6.t(k)   = t_sol(end);
    X.n7.y(k)   = y_sol(end,23);
    X.n7.t(k)   = t_sol(end);
    X.n8.y(k)   = y_sol(end,24);
    X.n8.t(k)   = t_sol(end);
    X.n9.y(k)   = y_sol(end,25);
    X.n9.t(k)   = t_sol(end);
    X.nm.y(k)   = y_sol(end,26);
    X.nm.t(k)   = t_sol(end);
    X.phi0.y(k) = y_sol(end,27);
    X.phi0.t(k) = t_sol(end);
    X.CO2outgas.y(k)   = y_sol(end,28);
    X.CO2outgas.t(k)   = t_sol(end);
    X.CO2_d.t(k) = t_sol(end);
    X.CO2_d.y(k) = y_sol(end,29);
    X.PAA.y(k)   = y_sol(end,30);
    X.PAA.t(k)   = t_sol(end);
    X.NH3.y(k)   = y_sol(end,31);
    X.NH3.t(k)   = t_sol(end);
    X.mu_P_calc.y(k)   = y_sol(end,32);
    X.mu_P_calc.t(k)   = t_sol(end);
    X.mu_X_calc.y(k)   = y_sol(end,33);
    X.mu_X_calc.t(k)   = t_sol(end);
%     X.S_pred.t(k)   = t_sol(end);
    X.X.y(k)   = X.a0.y(k) + X.a1.y(k) + X.a3.y(k) + X.a4.y(k);
    X.X.t(k) = t_sol(end);    
    X.Fault_ref.y(k) = u.Fault_ref; 
    X.Fault_ref.t(k) = t_sol(end);
    X.Control_ref.y(k) =  Ctrl_flags.PRBS;
    X.Control_ref.t(k) = Ctrl_flags.Batch_Num; 
    X.PAT_ref.y(k) = Ctrl_flags.Raman_spec;
    X.PAT_ref.t(k) = Ctrl_flags.Batch_Num; 
    X.Batch_ref.t(k) = Ctrl_flags.Batch_Num; 
    X.Batch_ref.y(k) = Ctrl_flags.Batch_Num; 
  
O2_in = 0.204; % % oxygen  in air
%% Calculating the OUR / CER
X.OUR.y(k) = (32*X.Fg.y(k)/22.4)*(O2_in-X.O2.y(k)*(0.7902/(1-X.O2.y(k)-X.CO2outgas.y(k)/100)));% in grams of oxygen consumed per sec
X.OUR.t(k) =  t_sol(end); 
% Calculating the CER
X.CER.y(k) = (44*X.Fg.y(k)/22.4)*((0.65*X.CO2outgas.y(k)/100)*(0.7902/(1-O2_in-X.CO2outgas.y(k)/100)-  0.0330));
X.CER.t(k) = t_sol(end);
% Adding in Raman Spectra
if k>10
if Ctrl_flags.Raman_spec==1;    
X =Raman_Sim(k,X,h, T);
elseif Ctrl_flags.Raman_spec==2
     X =Raman_Sim(k,X,h, T);
     X = Substrate_prediction(k,X,h,T);
end 
end
%% Off-line measurements recorded


if rem(t_sol(end),Ctrl_flags.Off_line_m)==0 || t_sol(end) ==1 || t_sol(end) ==T
    delay = Ctrl_flags.Off_line_delay;
    X.NH3_offline.y(k)   = X.NH3.y(k-delay);
    X.NH3_offline.t(k)   = X.NH3.t(k-delay);
    X.Viscosity_offline.y(k)   = X.Viscosity.y(k-delay);
    X.Viscosity_offline.t(k)   = X.Viscosity.t(k-delay);
    X.PAA_offline.y(k)   = X.PAA.y(k-delay);
    X.PAA_offline.t(k)   = X.PAA.t(k-delay);
    X.P_offline.y(k)   = X.P.y(k-delay);
    X.P_offline.t(k)   = X.P.t(k-delay);
    X.X_offline.y(k)   = X.X.y(k-delay)';
    X.X_offline.t(k)   = X.X.t(k-delay)';
else 
    X.NH3_offline.y(k)   = NaN;
    X.NH3_offline.t(k)   = NaN;
    X.Viscosity_offline.y(k)   = NaN;
    X.Viscosity_offline.t(k)   = NaN;
    X.PAA_offline.y(k)   = NaN;
    X.PAA_offline.t(k)   = NaN;
    X.P_offline.y(k)   = NaN;
    X.P_offline.t(k)   = NaN; 
    X.X_offline.y(k)   = NaN';
    X.X_offline.t(k)   = NaN';
end

  

end 
% close(h)
% unit conversions
X.pH.y = -log(X.pH.y) / log(10); % convert to pH from H+ concentration
X.Q.y = X.Q.y / 1000; % convert heat from Qrxn to kcal
%X.CO2.yUnit = '% outgas';
%% Add sensor noise 
% Add sensor noise to outputs after simulation can also be added throughout
% batch by including in ode loop and also can be applied to other variables
% T,pH,O2..
%  X.DO2.y = X.DO2.y + randn(N,1)*0.01;
%  X.CO2outgas.y = X.CO2outgas.y + randn(N,1)*0.01;
 end

%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037


