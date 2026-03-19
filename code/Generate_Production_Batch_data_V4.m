%% Generate multiple batches for a given Production phase. 
% IndPenSim_V2.01 Main file
% Last editted 9th of July 2019
% Objective: Generates batches for analysis
%% Copyright
% Stephen Goldrick Apr 2019 
% Contact: s.goldrick@ucl.ac.uk
% University College London, The University of Manchester, Newcastle University and Perceptive Engineering
% Please reference: "The Development of an Industrial Scale Fed-Batch
% Fermentation Simulation", Stepen Goldrick, Andrei Stefen, David Lovett,
% Gary Montague, Barry Lennox, Journal of Biotechnology 2015
%(https://www.sciencedirect.com/science/article/pii/S0168165614009377)
% And also please reference:  
% and "Modern day control challenges for industrial-scale fermentation
% processes" Goldrick et al. 2018 Computers and Chemical Engineering 
% All rights reserved. Copyright (c) The University of Manchester, University College London, Newcastle University and Perceptive Engineering.
% Clear workspace
 close all 
 clear all 
 clc

Data_generation_flag= 2; % 1-  Generate data for a fixed time interval
                         % 2- Generate Fixed number of batches (generate
                         % Normal and Fault batches
Operational_days = 336; 
Bioreactor_turn_around_time =3;
if  Data_generation_flag==1
% Select number of batches to generate
Production_Phase_in_years = 1; % Selecting a Production phase in year (1 year =  365 days)
Max_theoritical_number_of_batches =  round((Production_Phase_in_years*Operational_days)/11); 
Num_of_Batches = Max_theoritical_number_of_batches; % Number of Batches to generate
Batch_fault_order_reference= [0*ones(Num_of_Batches,1)]; 
else  
    Batch_run_flags.Batch_fault_order_reference = [0,1]; % Fault reference 
    Batch_run_flags.Control_strategy = [0,1]; % 0 - Recipe driven (i.e Sequential batch control (SBC)) 
                                            % 1- Operator conntroller batches   
    Batch_run_flags.Batch_length = [1,0]; % 0 - Fixed Batch length
                                          % 1 - Uneven batch length
    Batch_run_flags.Raman_spec =[1,2];  % 0 - Don't Record Raman data 
                                          % 1 - Record Raman Data
                                          % 2-  Use Raman data to control PAA
    Num_of_Batches = size(Batch_run_flags.Batch_fault_order_reference,2);
    Production_Phase_in_years =  (Num_of_Batches*(11+3))/Operational_days; 
    
end 
Save_batch_flag = 1; % 1- save data 0 - don't save data

Batches_File_Name  = 'IndPenSim_V2_export_V6'; % File name for saving generated batches
% Flags outlining Batch Operation
Generate_batch_flag = 1;% 1- Generate Batches and Plot figures ,  0 -  Load data from workspace and plot figures
% Batch operation adjustments to made in indpensim.m
if Generate_batch_flag == 1
Batch_start = 1;
for Batch_no = Batch_start:1:Num_of_Batches      
     matFileName = sprintf('Batch_%02d', Batch_no);  
     Xref = indpensim_run(Batch_no, Batch_run_flags); 
     Raw_Batch_data.(matFileName) = Xref;   
     Summmary_of_campaign(Batch_no,:)  = [Raw_Batch_data.(matFileName).Stats.Penicllin_harvested_during_batch , Raw_Batch_data.(matFileName).Stats.Penicllin_harvested_end_of_batch,Raw_Batch_data.(matFileName).Stats.Penicllin_yield_total,  Raw_Batch_data.(matFileName).Fg.t(end)/24]; 
     Summary_of_batch_lengths  = sum(ceil(Summmary_of_campaign(:,4))) +Batch_no*Bioreactor_turn_around_time; 
     if Summary_of_batch_lengths >Production_Phase_in_years*Operational_days
         break 
     else 
         
     end 
end
Num_of_Batches = Batch_no; 




%% 
Batch_Records= Generate_Batch_records(Raw_Batch_data,Batches_File_Name, Batch_run_flags); 
if Save_batch_flag ==1
savename = Batches_File_Name; 
save(savename);
clearvars -except Batch_Records Num_of_Batches Batch_start Summmary_of_campaign Raw_Batch_data
end
Generate_CSV_file_flag =1;
if Generate_CSV_file_flag ==1
    All_variables = 1; 
    
    
    
end 

end 


All_variables_names =fieldnames(Batch_Records.Batch_01);
Var_all = size(All_variables_names, 1);

IndPenSim_QbD_Figure_properties; 
set(0,'DefaultFigureWindowStyle','docked') 


for plot_no = 1:1:Var_all(end)
    % Defining the figure properties
figure
  hold on 
        set(gcf, 'color', [1 1 1])
        set(gca,'color', [1 1 1])
       
grid on
Temp_Var = All_variables_names(plot_no);
Temp_Var_name = regexprep(Temp_Var, '_',' ', 'emptymatch');

    
 Temp_title = (Temp_Var_name);

 set(gcf,'name',char(Temp_title))  
hold on
% Plotting all batches 
 for Batch_no = Batch_start:1:Num_of_Batches;        
    matFileName = sprintf('Batch_%02d', Batch_no); 
    Var_name_temp = char(Temp_Var);  
if  isfield(Batch_Records.(matFileName).(Var_name_temp), 't') ==1
    Temp_T_value  =  Batch_Records.(matFileName).(Var_name_temp).t;
    Temp_Y_value =   Batch_Records.(matFileName).(Var_name_temp).y;
   plot(Temp_T_value, Temp_Y_value,'Color',cmap(Batch_no,:),'lineStyle','-','linewidth',Line_width_fig);
     if diff(Temp_T_value) ==0; 
       plot(Temp_T_value(1), Temp_Y_value(1),'d', 'MarkerSize', 20,'MarkerFaceColor', 'r'   )
     
   end 
   xlabel(['Time [' Batch_Records.(matFileName).(Var_name_temp).tUnit ']'] ,'fontsize',Font_size_fig1 );
ylabel([Batch_Records.(matFileName).(Var_name_temp).name ' [' Batch_Records.(matFileName).(Var_name_temp).yUnit ']'] ,'fontsize',Font_size_fig1 );
title([Batch_Records.(matFileName).(Var_name_temp).name ' plot'],'fontsize',Font_size_fig1 );
   
else 
    close(gcf)    
end 
 end
end

% % Plot spectral data of final batch
if ismember('Raman_Spec',All_variables_names)==1
figure
 set(gcf,'name','Raman Spectra')  
  hold on 
        set(gcf, 'color', [1 1 1])
        set(gca,'color', [1 1 1])    
grid on
plot(Batch_Records.(matFileName).Raman_Spec.Wavelength, Batch_Records.(matFileName).Raman_Spec.Intensity)
 set(get(gca,'XLabel'),'String','Wavelength (cm^{-1}','fontsize',Font_size_fig1 )
 set(get(gca,'YLabel'),'String','Intensity (a.u)','fontsize',Font_size_fig1)
 set(get(gca,'title'), 'String','Raman spectra','fontsize',Font_size_fig1)
end



%% Copyright
% Stephen Goldrick 
% University of Manchester, Newcastle University and Perceptive Engineering and University College London
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

% Please reference: DOI: 10.1016/j.jbiotec.2014.10.029 and https://doi.org/10.1016/j.compchemeng.2019.05.037




