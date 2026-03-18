%% IndPenSim_QbD  Figure properties
% Defining Figure properties for Gluc vs OTR paper 
% Defining figure properties#
IndPenSim_QbD_paper =1; 
if IndPenSim_QbD_paper ==1
axis_pos = 0.95; 
Font_size_fig1 = 14; 
Font_size_fig2 = 10;
cmap = distinguishable_colors(Num_of_Batches);
Color_reference = [1,2,3,4,5,6,7,8,9,10];
plotStyle = {'-.',':','-','--','-','-.'}; % add as many as you need
plot_LS = plotStyle(1);
MarkerStyle = {'none', 'none', 'none', 'none','>','o' };
Temp_Marker_size =6;
Line_width_fig = 2; 
lengend_font_size= 6;
Marker_spread = 40; 
 %set(0,'DefaultFigureWindowStyle','docked') 
 set(0,'DefaultFigureWindowStyle','normal')
 
   MVDA_legend_font   =28;
MVDA_Font_size_fig = 28;
MVDA_line_width_border = 3;
MVDA_line_width_fig = 3;
MVDA_export_paper =0   ;
MVDA_line_width_gscatter = 2;
MVDA_Font_size_fig_gscatter=12;
 
end 

