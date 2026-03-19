%Generate Batch records from Output data
function Batch_Records= Generate_Batch_records(Raw_Batch_data,Batches_File_Name, Batch_run_flags )
%Batch_Records= Generate_Batch_records(Raw_Batch_data)

All_batches = fieldnames(Raw_Batch_data); 
Num_of_Batches = size(All_batches,1); 
Temp_batch = Raw_Batch_data.(All_batches{1}); 
%Temp_batch  = Raw_Batch_data.Batch_01; 
%field = {'sc', 'abc', 'a0','a1','a3','a4','n0','n1', 'n2','n3','n4','n5','n6','n7','n8','n9','nm','phi0','Culture_age','mup', 'mux', 'X_CER', 'mu_X_calc','mu_P_calc','F_discharge_cal', 'CO2_d','S_pred', 'NH3', 'PAA', 'Viscosity', 'X'};

field = {'sc', 'abc', 'a0','a1','a3','a4','n0','n1', 'n2','n3','n4','n5','n6','n7','n8','n9','nm','phi0','Culture_age','mup', 'mux', 'X_CER', 'mu_X_calc','mu_P_calc','F_discharge_cal', 'CO2_d', 'NH3', 'PAA', 'Viscosity', 'X','PRBS_noise_addition' };
%  Temp_batch = rmfield(Raw_Batch_data.(Batch_ref),field);


 % Generate csv file
Temp_batch = rmfield(Raw_Batch_data.Batch_01,field);
All_variables_to_import = fieldnames(Temp_batch); 

    

%if ismember('Raman_Spec',All_variables_to_import)==1

All_Batch_stats = [];
  DF_All_Batch_data = [];
 for  Batch_no = 1:1:Num_of_Batches
     Batch_ref = sprintf('Batch_%02d', Batch_no); 
     Temp_time = Raw_Batch_data.(Batch_ref).(All_variables_to_import{1}).t; 
     DF_Batch_data = Temp_time; 
     DF_all_data_headers =  {Raw_Batch_data.(Batch_ref).(All_variables_to_import{1}).name}; 
%      All_variables_to_import = fieldnames(Raw_Batch_data.(Batch_ref)); 
     if ismember('PAA_pred',All_variables_to_import)
         field = {'sc', 'abc', 'a0','a1','a3','a4','n0','n1', 'n2','n3','n4','n5','n6','n7','n8','n9','nm','phi0','Culture_age','mup', 'mux', 'X_CER', 'mu_X_calc','mu_P_calc','F_discharge_cal', 'CO2_d', 'NH3', 'PAA', 'Viscosity', 'X','PRBS_noise_addition', 'PAA_pred' };
   Temp_batch = rmfield(Raw_Batch_data.(Batch_ref),field);
    All_variables_to_import = fieldnames(Temp_batch); 
    end 
     
     for ii = 1:1:(size(All_variables_to_import,1) -2)
     DF_Batch_data  =  [DF_Batch_data, Raw_Batch_data.(Batch_ref).(All_variables_to_import{ii}).y];           
     end 
     if sum(Raw_Batch_data.(Batch_ref).Fault_ref.y) ==0
         Batch_fault_ref = 0*ones(size(DF_Batch_data,1),1);
     else
         Batch_fault_ref = 1*ones(size(DF_Batch_data,1),1);
     end
     Temp_stat_data = [-Raw_Batch_data.(Batch_ref).Stats.Penicllin_harvested_during_batch/1000, Raw_Batch_data.(Batch_ref).Stats.Penicllin_harvested_end_of_batch/1000, Raw_Batch_data.(Batch_ref).Stats.Penicllin_yield_total/1000];
     Raman_data = [];
      if ismember('Raman_Spec',All_variables_to_import)==1
    DF_Batch_data = [DF_Batch_data, ones(size(DF_Batch_data,1),1)*Batch_no,Batch_fault_ref , Raw_Batch_data.(Batch_ref).Raman_Spec.Intensity'];
      else
          DF_Batch_data = [DF_Batch_data, Batch_fault_ref,  ones(size(DF_Batch_data,1),1)*Batch_no   ];
      end 
    DF_All_Batch_data = [DF_All_Batch_data;DF_Batch_data] ;  
     All_Batch_stats = [All_Batch_stats; Batch_no, Temp_stat_data, Batch_fault_ref(end) ]
 end
 
 All_variables_to_import(1:end-2)
 
 for oo =1:1:size(All_variables_to_import(1:end-2),1)
     Temp_variable = strcat(Raw_Batch_data.Batch_01.(All_variables_to_import{oo}).name, '(',  All_variables_to_import{oo}, ':', Raw_Batch_data.Batch_01.(All_variables_to_import{oo}).yUnit,')')
 All_variables_with_units{oo} = Temp_variable;
 end 
 
 
 if ismember('Raman_Spec',All_variables_to_import)==1 
Raman_wavelength_headers = string(Raw_Batch_data.(Batch_ref).Raman_Spec.Wavelength)';
df_headers = cellstr([{'Time (h)'}; All_variables_with_units'; {'Batch ID'}; {'Fault flag'};  Raman_wavelength_headers']);
 else
     df_headers = [{'Time (h)'}; All_variables_to_import(1:end-1); {'Batch_ref'};  ];
 end 

CSV_file_temp = strcat(Batches_File_Name, '.csv')
csvwrite_with_headers(CSV_file_temp, DF_All_Batch_data,df_headers')




df_Stats_headers = [{'Batch ref'; 'Penicllin_harvested_during_batch(kg)' ; 'Penicllin_harvested_end_of_batch (kg)'; 'Penicllin_yield_total (kg)'; 'Fault ref(0-NoFault 1-Fault)' }];
CSV_file_temp_stats = strcat(Batches_File_Name, '_Statistics','.csv')
csvwrite_with_headers(CSV_file_temp_stats, All_Batch_stats,df_Stats_headers )
 
 


 for Batch_no = 1:1:Num_of_Batches;        
     matFileName = sprintf('Batch_%02d', Batch_no); 
     [row] = find(~isnan(Raw_Batch_data.(matFileName).PAA_offline.y));
     Raw_Batch_data.(matFileName).PAA_offline.y = Raw_Batch_data.(matFileName).PAA_offline.y(row);
     Raw_Batch_data.(matFileName).PAA_offline.t = Raw_Batch_data.(matFileName).PAA_offline.t(row);
     Raw_Batch_data.(matFileName).P_offline.y = Raw_Batch_data.(matFileName).P_offline.y(row);
     Raw_Batch_data.(matFileName).P_offline.t = Raw_Batch_data.(matFileName).P_offline.t(row);
     Raw_Batch_data.(matFileName).NH3_offline.y = Raw_Batch_data.(matFileName).NH3_offline.y(row);
     Raw_Batch_data.(matFileName).NH3_offline.t = Raw_Batch_data.(matFileName).NH3_offline.t(row);  
     Raw_Batch_data.(matFileName).X_offline.y = Raw_Batch_data.(matFileName).X_offline.y(row)';
     Raw_Batch_data.(matFileName).X_offline.t = Raw_Batch_data.(matFileName).X_offline.t(row)';   
     Raw_Batch_data.(matFileName).Viscosity_offline.y = Raw_Batch_data.(matFileName).Viscosity_offline.y(row);
     Raw_Batch_data.(matFileName).Viscosity_offline.t = Raw_Batch_data.(matFileName).Viscosity_offline.t(row);   
     Raw_Batch_data.(matFileName).O2.y = Raw_Batch_data.(matFileName).O2.y*100; % Convert to %  
%      Raw_Batch_data.(matFileName).Foil.t = Raw_Batch_data.(matFileName).Foil.t';   
%      Raw_Batch_data.(matFileName).Foil.y = Raw_Batch_data.(matFileName).Foil.y';  
%         Raw_Batch_data.(matFileName).Fpaa.t = Raw_Batch_data.(matFileName).Fpaa.t';   
%        Raw_Batch_data.(matFileName).Fpaa.y = Raw_Batch_data.(matFileName).Fpaa.y';     
     
     Temp_batch = rmfield(Raw_Batch_data.(matFileName),field) ;
     Batch_Records.(matFileName) = Temp_batch;          
 end 
  
  

 save_file_name = Batches_File_Name; 
 save(save_file_name, 'Batch_Records');
