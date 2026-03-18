% Generating Simulated Raman Spectra
function    X = Raman_Sim(k,X,h, T) 

   if (rem(k,1)==0) % Build Spectra   
% Building history of Raman Spectra 
    Wavenumber_max = 2200;
    Intensity_shift = ones([Wavenumber_max,1]);
Intensity_shift1= ones([Wavenumber_max,1]);
for j = 1:1:Wavenumber_max
     b = j/(Wavenumber_max*0.5);
     Intensity_shift1(j,1) = exp(b) -0.5;
 
end 
Intensity_shift = ones([Wavenumber_max,1]);
 New_Spectra = ones([Wavenumber_max,1]);


a = -0.00178143846614472*0.1;
b = 1.05644816081515;
c = -0.0681439987249108*0.1; 
d = -0.02;
Product_S = X.P.y(k)/40; %Penicillin Scaled 
Biomass_S =X.X.y(k)/40; % Biomass Scaled
Viscosity_S = X.Viscosity.y(k)/100; 
Time_S = k/(T/h);
Intensity_increase1(k) = a*Biomass_S+ b*Product_S+ c*Viscosity_S +d*Time_S;% Spectra_no = 1:1:5; 
scaling_factor =   370000; 
Gluc_increase = 800000*3/1400;

PAA_increase=  1700000/1000;
Prod_increase = 100000;

  %% Loading in the reference Raman Spectral file

  reference_Spectra_2200 = load('reference_Specra.txt');

 
   X.Raman_Spec.Wavelength  = reference_Spectra_2200(1:Wavenumber_max,1); 
    reference_spectra = reference_Spectra_2200(1:Wavenumber_max,2);   
      New_Spectra(:,k) = Intensity_increase1(k).*Intensity_shift1.*scaling_factor +reference_spectra;  
X.Raman_Spec.Intensity(:,k)  = New_Spectra(:,k);
    
random_noise = ones([(Wavenumber_max+1),1]);
one_matrix = ones([Wavenumber_max+1,1]);
random_noise_summed = ones([Wavenumber_max,1]);
New_Spectra_noise = ones([Wavenumber_max,571]);

random_number = randi(3,Wavenumber_max,1);


for i = 1:1:(Wavenumber_max)
    noise_factor =50; 
if random_number(i) == 1; 
    random_noise(i) = 0;
elseif random_number(i) == 2;
    random_noise(i) = noise_factor;
else 
    random_noise(i) = -noise_factor;
end 
end 
for i = 1:1:Wavenumber_max
random_noise_summed(i) = sum(random_noise(1:i)); 
end 
random_noise_summed_smooth = smooth(random_noise_summed,25);




 New_Spectra_noise(:,k)= New_Spectra(:,k) +10*random_noise_summed_smooth;
 
 X.Raman_Spec.Intensity(:,k)  = New_Spectra_noise(:,k);
 


 % -----------------------------------------------------------------------
%             Aim 3. Creating the bell curve response for Glucose
%----------------------------------------------------------------------

Glucose_raw_peaks = ones([Wavenumber_max,1])*0;
Glucose_raw_peaks_G_peaka = ones([Wavenumber_max,1])*0;
Glucose_raw_peaks_G_peakb = ones([Wavenumber_max,1])*0;
Glucose_raw_peaks_G_peakc = ones([Wavenumber_max,1])*0;
PAA_raw_peaks_G_peaka = ones([Wavenumber_max,1])*0;
PAA_raw_peaks_G_peakb = ones([Wavenumber_max,1])*0;
Product_raw_peaka = ones([Wavenumber_max,1])*0;
Product_raw_peakb = ones([Wavenumber_max,1])*0;

%% Glucose peaks
%---------Peak A  
 peaka = 219;
 peaka_width = 70;
     peaka_lenght  = peaka_width*2; 
     peaka_std_dev = peaka_width/2;
     mean = 0; 
     for x = -peaka_lenght:1:peaka_lenght
     Glucose_raw_peaks_G_peaka(x+peaka) = (peaka_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peaka_std_dev).^2);
     end 
%---------- Peak B 

  peakb = 639;
 peakb_width = 20;
     peakb_lenght  = peakb_width*2; 
     peakb_std_dev = peakb_width/2;
     mean = 0; 
     for x = -peakb_lenght:1:peakb_lenght
     Glucose_raw_peaks_G_peakb(x+peakb) = ((peakb_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peakb_std_dev).^2)/4.3);
     end 

%------------ Peak C
  peakc = 1053;
 peakc_width = 100;
     peakc_lenght  = peakc_width*2; 
     peakc_std_dev = peakc_width/2;
     mean = 0; 
     for x = -peakc_lenght:1:peakc_lenght
     Glucose_raw_peaks_G_peakc(x+peakc) = (peakc_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peakc_std_dev).^2);
     end
     
%% PAA  peaks
%---------Peak A  
 peaka = 419;
 peaka_width = 60;
     peaka_lenght  = peaka_width*2; 
     peaka_std_dev = peaka_width/2;
     mean = 0; 
     for x = -peaka_lenght:1:peaka_lenght
     PAA_raw_peaks_G_peaka(x+peaka) = (peaka_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peaka_std_dev).^2);
     end 
%---------- Peak B 

  peakb = 839;
 peakb_width = 15;
     peakb_lenght  = peakb_width*2; 
     peakb_std_dev = peakb_width/2;
     mean = 0; 
     for x = -peakb_lenght:1:peakb_lenght
     PAA_raw_peaks_G_peakb(x+peakb) = ((peakb_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peakb_std_dev).^2)/4.3);
     end 
     
% Adding in  Peak aPen G Peak
     
      peakPa = 800;
 peakPa_width = 30;
     peakPa_lenght  = peakPa_width*4; 
     peakPa_std_dev = peakPa_width/2;
     mean = 0; 
     for x = -peakPa_lenght:1:peakPa_lenght
     Product_raw_peaka(x+peakPa) = (peakPa_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peakPa_std_dev).^2);
     end 
       % Adding in  Peak b for Pen G Peak
     
      peakPb = 1200;
     peakPb_width = 30;
     peakPb_lenght  = peakPb_width*30; 
     peakPb_std_dev = peakPb_width/2;
     mean = 0; 
     for x = -peakPb_lenght:1:peakPb_lenght
     Product_raw_peakb(x+peakPb) = (peakPb_std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/peakPb_std_dev).^2);
     end 
   
 
     
     total_peaks_G = Glucose_raw_peaks_G_peaka+Glucose_raw_peaks_G_peakb+Glucose_raw_peaks_G_peakc;
     total_peaks_PAA = PAA_raw_peaks_G_peaka+PAA_raw_peaks_G_peakb; 
     total_peaks_P = Product_raw_peakb + Product_raw_peaka; 
 K_G =0.005;
 K_PAA = 4000;
 Substrate_raman = X.S.y(k);
 PAA_raman = X.PAA.y(k);
 
  New_Spectra_noise_Gluc_Prod(:,k) = New_Spectra_noise(:,k)+ total_peaks_G*Gluc_increase* Substrate_raman/(K_G+ Substrate_raman) +total_peaks_PAA*PAA_increase*PAA_raman +total_peaks_P*Prod_increase.*X.P.y(k);    

 X.Raman_Spec.Intensity(:,k)  = New_Spectra_noise_Gluc_Prod(:,k);
 
   end 
