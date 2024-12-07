function [FILT]= filtering (fsamp,sd_sig,muscle_code)

% Filtro passabanda per la rimozione del rumore
% High Pass Filter (LPF)
[n_chs n_samps]= size(sd_sig.GM);
fNy=fsamp/2;

Wp_high = 20/fNy;  
Ws_high = 10/fNy;   
Rp_high = 5;  % dB
Rs_high = 40;  % dB

[Order_high,Wn_high] = buttord(Wp_high,Ws_high,Rp_high,Rs_high); 
[b_high,a_high] = butter(Order_high,Wn_high,"high"); 

% Low Pass Filter (LPF):
Wp_low = 400/fNy;  
Ws_low = 450/fNy;   
Rp_low = 5;  % dB
Rs_low = 20;  % dB

[Order_low,Wn_low] = buttord(Wp_low,Ws_low,Rp_low,Rs_low);
[b,a] = butter(Order_low,Wn_low,"low");              

% Filtri per l'estrapolazione dell'inviluppo
% High Pass Filter (HPF)
Wp_h = 35/fNy;  
Ws_h = 10/fNy;   
Rp_h = 5;  % dB
Rs_h = 20;  % dB

[Order_h,Wn_h] = buttord(Wp_high,Ws_high,Rp_high,Rs_high); 
[b_h,a_h] = butter(Order_h,Wn_h,"high"); 

% Low Pass Filter (HPF)  
Order_lowi=4;
Wn_lowi1=3.8/fNy;
Wn_lowi2=9/fNy;
[bi1,ai1] = butter(Order_lowi,Wn_lowi1,"low");    
[bi2,ai2] = butter(Order_lowi,Wn_lowi2,"low");   

FILT=struct;
for i=1:length(muscle_code)
  FILT.(muscle_code{i})=[];
end
for i=1:length(muscle_code)
    if i<length(muscle_code)
         
        if i<=3
         % RIMOZIONE RUMORE
         % HPF Butterworth
         tmp=sd_sig.(muscle_code{i})-mean(sd_sig.(muscle_code{i}),'omitnan'); % rimozione valor medio
         tmp=fillmissing(tmp,'pchip');
         tmpEMGfilt = filtfilt(b_high,a_high,tmp); 
         % LPF Butterworth
         tmpEMGfilt_final = filtfilt(b,a,tmpEMGfilt);  

         % INVILUPPO
         % HPF Butterworth
         EMGfilt_final=filtfilt(b_h,a_h,tmpEMGfilt_final);
         % Full-wave rectification
         EMGfiltrect = abs(EMGfilt_final);
         % LPF Butterworth
         EMGfiltrect_final=filtfilt(bi1,ai1,EMGfiltrect);
         FILT.(muscle_code{i})(1,:)=EMGfiltrect_final;  
        else
         % RIMOZIONE RUMORE
         % HPF Butterworth
         tmp=sd_sig.(muscle_code{i})-mean(sd_sig.(muscle_code{i}),'omitnan'); % rimozione valor medio
         tmp=fillmissing(tmp,'pchip');
         tmpEMGfilt = filtfilt(b_high,a_high,tmp); 
         % LPF Butterworth
         tmpEMGfilt_final = filtfilt(b,a,tmpEMGfilt);  

         % INVILUPPO
         % HPF Butterworth
         EMGfilt_final=filtfilt(b_h,a_h,tmpEMGfilt_final);
         % Full-wave rectification
         EMGfiltrect = abs(EMGfilt_final);
         % LPF Butterworth
         EMGfiltrect_final=filtfilt(bi2,ai2,EMGfiltrect);
         FILT.(muscle_code{i})(1,:)=EMGfiltrect_final;
        end
            
    else
        for i_ch=1:n_chs
         % RIMOZIONE RUMORE
         % HPF Butterworth
         tmp=sd_sig.(muscle_code{i})(i_ch,:)-mean(sd_sig.(muscle_code{i})(i_ch,:)); %rimozione valor medio
         tmp=fillmissing(tmp,'pchip'); %interpolazione per rimuovere i NaN
         tmpEMGfilt = filtfilt(b_high,a_high,tmp); 
        
         % LPF Butterworth:
         tmpEMGfilt_final = filtfilt(b,a,tmpEMGfilt);  

         % INVLUPPO
         % HPF Butterworth
         EMGfilt_final=filtfilt(b_h,a_h,tmpEMGfilt_final);
         % Full-wave rectification
         EMGfiltrect = abs(EMGfilt_final);
         % LPF Butterworth
         EMGfiltrect_final=filtfilt(bi2,ai2,EMGfiltrect);
         FILT.(muscle_code{i})(i_ch,:)=EMGfiltrect_final;   
        end
      
    end
        tmp=[];
        tmpEMGfilt=[];
        EMGfiltrect=[];
        EMGfiltrect_final =[];
end

end
