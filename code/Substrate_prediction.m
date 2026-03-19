function X = Substrate_prediction(k,X,h,T)

load PAA_PLS_model
 j = k-1;
 
 options1 = [];
Raman_Spec_sg = sgolayfilt(X.Raman_Spec.Intensity(:,j)',2,5);
Raman_Spec_sg  = Raman_Spec_sg';
Raman_Spec_sg_d = diff(Raman_Spec_sg);
PAA_peaks_Spec = Raman_Spec_sg_d([350:500 800:860],:);
 No_LV =4;
 X.PAA_pred.y(j) = PAA_peaks_Spec'*b(No_LV,:)';
 if j > 20 
     X.PAA_pred.y(j) = (X.PAA_pred.y(j-1) + X.PAA_pred.y(j-2) +X.PAA_pred.y(j))/3;
 end 
 

