clear all
close all
clc
% Studenti:
% Mazzoletti Laura 292229
% Richetto Francesco 292281
% Sabbadini Benedetta 292182
% Scimeca Sabrina 279968
% Scorrano Luca 296047

%% define global variables to reduce parameters in the functions
global markers
global joints
global order_pelvis
global left
global right
global muscle_code
global angles
global ref

%% Loading Data
kin_fsamp=100; %sample frequency

emg_fsamp= 2048; %sampling frequency
muscle_code=["BF","GMed","TA","VL","SOLm", "GM"];
markers= {'RASI'; 'LASI'; 'RPSI'; 'LPSI';'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'; ...
    'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'};
joints={'PC'; 'RK'; 'LK'; 'RA'; 'LA'; 'RTOE';'LTOE'};
angles={'hip'; 'knee'; 'ankle';'hip_r'; 'knee_r'; 'ankle_r'};
order_pelvis={'LASI';'RASI';'RPSI';'LPSI';'LASI'};
left={'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'};
right={'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'};
ref={'pelvis','femur_r','femur_l','thigh_r','thigh_r_nt','thigh_l', 'thigh_l_nt','foot_l','foot_r'};


%Load marker coordinates
 [sd_sig, traj]= load_kin_EMG_data();

%swap columns 1 and 2 in order to have medio lateral axis along y
for i_m=1:length(markers)
    traj.(markers{i_m})= traj.(markers{i_m})(:,[2 1 3]);
end

% calcolo sistemi di riferimento
JOINTS=struct;
loc_ref=struct;
[loc_ref.pelvis,JOINTS.PC,JOINTS.RH,JOINTS.LH]=calc_references_pelvis(traj);
[loc_ref.femur_r,loc_ref.femur_l,JOINTS.RK,JOINTS.LK]=calc_references_femur(traj,JOINTS);
[loc_ref.thigh_r, loc_ref.thigh_r_nt,loc_ref.thigh_l, loc_ref.thigh_l_nt,JOINTS.RA, JOINTS.LA]= calc_references_shank(traj,JOINTS,loc_ref.femur_r,loc_ref.femur_l);
[loc_ref.foot_r,loc_ref.foot_l,JOINTS.RTOE ,JOINTS.LTOE]= calc_references_foot(traj,JOINTS);



%% calcolo angoli articolari

ANGLES=struct;
% LEFT
for i=1:size(loc_ref.pelvis,1)
[ANGLES.hip(i,:)]=joint_angles([loc_ref.pelvis(i,1:3)' loc_ref.pelvis(i,7:9)' loc_ref.pelvis(i,4:6)'] ,[loc_ref.femur_l(i,1:3)' loc_ref.femur_l(i,7:9)' loc_ref.femur_l(i,4:6)']);
[ANGLES.knee(i,:)]=joint_angles(-[loc_ref.femur_l(i,1:3)' loc_ref.femur_l(i,7:9)' loc_ref.femur_l(i,4:6)'],[loc_ref.thigh_l_nt(i,1:3)' loc_ref.thigh_l_nt(i,7:9)' loc_ref.thigh_l_nt(i,4:6)']);
[ANGLES.ankle(i,:)]=joint_angles([loc_ref.thigh_l(i,7:9)' loc_ref.thigh_l(i,1:3)' loc_ref.thigh_l(i,4:6)'],[loc_ref.foot_l(i,1:3)' loc_ref.foot_l(i,7:9)' loc_ref.foot_l(i,4:6)']);
end
% % RIGHT
% for i=1:size(loc_ref.pelvis,1)
% [ANGLES.hip(i,:)]=joint_angles([loc_ref.pelvis(i,1:3)' loc_ref.pelvis(i,7:9)' loc_ref.pelvis(i,4:6)'] ,[loc_ref.femur_r(i,1:3)' loc_ref.femur_r(i,7:9)' loc_ref.femur_r(i,4:6)']);
% [ANGLES.knee(i,:)]=joint_angles([loc_ref.femur_r(i,1:3)' loc_ref.femur_r(i,7:9)' loc_ref.femur_r(i,4:6)'],[loc_ref.thigh_r_nt(i,1:3)' loc_ref.thigh_r_nt(i,7:9)' loc_ref.thigh_r_nt(i,4:6)']);
% [ANGLES.ankle(i,:)]=joint_angles(-[loc_ref.thigh_r(i,1:3)' loc_ref.thigh_r(i,7:9)' loc_ref.thigh_r(i,4:6)'],[loc_ref.foot_r(i,7:9)' loc_ref.foot_r(i,1:3)' loc_ref.foot_r(i,4:6)']);
% end




%% calcolo cicli del cammino
heel_l=traj.LHEE(:,3);
t=1/kin_fsamp:1/kin_fsamp:length(heel_l)/1/kin_fsamp;
maxi=max(heel_l);
[peaks,index]=findpeaks(-heel_l+maxi,'MinPeakDistance',108); 
index(1,:)=[];
peaks(1,:)=[];
% calcolo la distanza media dei passi per comprendere quando avviene la
% svolta
for i=1:length(index)-1
    dist(i)=index(i+1)-index(i);
end
meandist=mean(dist);

% prima colonna index:campione corrispondente al marker LHEE
% seconda colonna di index: maschera per eliminare il primo passo dopo ogni svolta

for i=1:length(index)-1
     if(index(i+1)-index(i))<meandist
         index(i+1,2)=0;
     else
       index(i+1,2)=1;
     end
end

%% Ricampionamento

samp= 100; 
for k=1:length(angles)/2
  ANGLES.(angles{k+3})=[];
    res_angle=[]; 
    for i=1:length(index)-1
        tmp=[];
        res_tmp=[];
        tmp=ANGLES.(angles{k})(index(i):index(i+1),1:3);
        res_tmp=resample(tmp,samp,size(tmp,1))';
        res_angle=[res_angle res_tmp];
       end
    ANGLES.(angles{k+3})(:,:)=res_angle';
end
steps=size(ANGLES.hip_r,1)/100;
%% esclusione dei passi non completi

TREND_ANGLES=struct;
nosteps=zeros(length(index),1);
% nosteps maschera per eliminare i passi dove gli angoli articolari hanno
% un NAN
for i=4:length(angles)
    for j=1:length(index)-1
        nan=0;
        for k=(j-1)*samp+1:j*samp
           nan=isnan(ANGLES.(angles{i})(k,:));
           for x=1:3
              if nan(x)==1
                 nosteps(j)=1;                
                break
              end
           end
        end
    end
end
% aggiungo la maschera dei nan alla maschera dei primi passi così da avere
% una maschera completa per l'eliminazione  
for i=1:length(nosteps)
    if (nosteps(i)==1)
        for j=1:size(index,1)
            if (i==j)
                index(j+1,2)=1;
            end
        end
    end
end
 % duration é la matrice in cui la prima colonna indica il campione d'inizio passo, la
 % seconda il campione di fine passo , la terza la durata in secondi e la
 % quarta l'indice corrispondente alla matrice index
duration=[];
k=0;
 for i=1:length(index)-1
     if(index(i,2)==0 && index(i+1,2)==0)
       k=k+1;
       duration(k,1)=index(i);
       duration(k,2)=index(i+1);
       duration(k,3)=(index(i+1)-index(i))/kin_fsamp;
       duration(k,4)=i;
      end
end
k=0;
for i=4:length(angles)
   TREND_ANGLES.(angles{i-3})=[];
   k=0;
   for j=1:size(duration,1)
      tmp=[];
      k=k+1;
      tmp=ANGLES.(angles{i})((duration(j,4)-1)*samp+1:duration(j,4)*samp,:);
      TREND_ANGLES.(angles{i-3})((k-1)*samp+1:k*samp,:)=tmp;
     
   end
end

steps=size(TREND_ANGLES.ankle,1)/samp;

figure()

plot(t,heel_l,t(duration(:,1)),heel_l(duration(:,1)),'bo'),title (' Cicli di cammino')
hold on
plot(t(duration(:,2)),heel_l(duration(:,2)),'r*'),title (' Cicli di cammino')
ylabel('Asse Z [mm]'), xlabel('Time [s]')

legend('','inizio passo','fine passo')

sprintf('Sono stati considerati %d passi',steps)
sprintf('La durata media di ogni passo è di %0.2f secondi',mean(duration(:,3)))

%% plot angoli articolari

close all

for i=1:length(angles)/2
    figure(i)
    abd_add=[];
    for j=2:steps 
        subplot(3,1,1)
        tmp=[];
        tmp=TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,1);
        abd_add=[abd_add tmp];
        plot(TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,1)),title("Abduzione-Adduzione"), axis tight
        xlabel('% of Gait'), ylabel('Gradi (°)')
        hold on
    end
    hold on

    intra_ext=[];
    for j=1:steps        
        subplot(3,1,2)  
        tmp=[];
        tmp=TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,2);
        intra_ext=[intra_ext tmp];
        plot(TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,2)),title("Intra-Extra Rotazione"), axis tight
        xlabel('% of Gait'), ylabel('Gradi (°)')
        hold on
    end
    hold on

    flex_ext=[];
    for j=1:steps
        subplot(3,1,3)
        tmp=[];
        tmp=TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,3);
        flex_ext=[flex_ext tmp];
        plot(TREND_ANGLES.(angles{i})((j-1)*samp+1:j*samp,3)),title("Flesso-Estensione"), axis tight
        xlabel('% of Gait'), ylabel('Gradi (°)')
        hold on
    end
    hold off
    sgtitle(["Angolo articolare: ",angles{i}]) 

    % plot intra_extra mediati con errore standard
    figure(4)
    subplot(3,3,(i*3)-2)
    mean_i_e=mean(intra_ext');
    std_i_e=std(intra_ext')/sqrt(length(mean_i_e));
    plot(0:100/length(mean_i_e+std_i_e):100*(length(mean_i_e+std_i_e)-1)/length(mean_i_e+std_i_e),mean_i_e+std_i_e,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_i_e-std_i_e):100*(length(mean_i_e-std_i_e)-1)/length(mean_i_e-std_i_e),mean_i_e-std_i_e,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_i_e):100*(length(mean_i_e)-1)/length(mean_i_e),mean_i_e,'b','LineWidth',1), grid on, title("Intra-Extra Rotazione : ",angles{i})
    xlabel('% of Gait'), ylabel('Gradi (°)')
    hold on 

    % plot abd-add mediati con errore standard
    subplot(3,3,(i*3)-1)
    mean_a_a=mean(abd_add');
    std_a_a=std(abd_add')/sqrt(length(mean_a_a));
    plot(0:100/length(mean_a_a+std_a_a):100*(length(mean_a_a+std_a_a)-1)/length(mean_a_a+std_a_a),mean_a_a+std_a_a,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_a_a-std_a_a):100*(length(mean_a_a-std_a_a)-1)/length(mean_a_a-std_a_a),mean_a_a-std_a_a,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_a_a):100*(length(mean_a_a)-1)/length(mean_a_a),mean_a_a,'b','LineWidth',1), grid on, title("Abduzione-Adduzione: ",angles{i})
    xlabel('% of Gait'), ylabel('Gradi (°)') 
    hold on   

    % plot flex-ext mediati con deviazione std
    subplot(3,3,i*3)
    mean_f_e=mean(flex_ext');
    std_f_e=std(flex_ext')/sqrt(length(mean_f_e));
    plot(0:100/length(mean_f_e+std_f_e):100*(length(mean_f_e+std_f_e)-1)/length(mean_f_e+std_f_e),mean_f_e+std_f_e,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_f_e-std_f_e):100*(length(mean_f_e-std_f_e)-1)/length(mean_f_e-std_f_e),mean_f_e-std_f_e,'r','LineWidth',0.5)
    hold on
    plot(0:100/length(mean_f_e):100*(length(mean_f_e)-1)/length(mean_f_e),mean_f_e,'b','LineWidth',1), grid on, title("Flesso-Estensione: ",angles{i})
    xlabel('% of Gait'), ylabel('Gradi (°)')
    hold off
   
end
 legend("errore standard sup","errore standard inf","angolo",'EdgeColor',[0 0 1]);

%% Filtraggio segnali muscolari

[FILT]= filtering (emg_fsamp,sd_sig,muscle_code);
[n_chs n_samps]= size(sd_sig.GM);
t= [1:n_samps]/emg_fsamp;

%% Plot segnali filtrati (inviluppo)
close all
vnorm_in1=[];
vnorm_in=[];
vnorm_in=max(max(abs(FILT.GM)));
figure()
for i=1:length(muscle_code)
    if i<length(muscle_code)
        subplot(5,1,i)
        vnorm_in1(i)=max(max(abs(FILT.(muscle_code{i})(1,:))));
        FILT.(muscle_code{i})(1,:)=abs(FILT.(muscle_code{i})(1,:));
        plot(t,FILT.(muscle_code{i})(1,:)/vnorm_in1(i),'k'), title (muscle_code{i});
        hold off
        sgtitle ('Inviluppi normalizzati'), xlabel('tempo (s)'),ylabel('EMG norm')
    else
        figure()
        for i_ch=1:n_chs
            FILT.(muscle_code{i})(i_ch,:)=abs(FILT.(muscle_code{i})(i_ch,:));
            plot(t,FILT.(muscle_code{i})(i_ch,:)/vnorm_in+n_chs-(n_chs-i_ch),'k'), title ('Inviluppi GM normalizzati'), xlabel('tempo (s)'), ylabel('EMG norm');
            hold on
        end
    end
end

%% calcolo ciclo passo per le attivazioni muscolari
duration=[];
k=0;
 for i=1:length(nosteps)-1
     if(nosteps(i)==0)
       k=k+1;
       duration(k,1)=index(i);
       duration(k,2)=index(i+1);
       duration(k,3)=index(i+1)-index(i);%dist(i);
      end
end
k=0;
for i=1:size(duration,1)
    if nosteps(i)==0
        k=k+1;
        start_i(k)=duration(i,1)/kin_fsamp; %trovo i tempi di inizio e fine passi 
        end_i(k)=duration(i,2)/kin_fsamp;
    end
end
for i=1:length(start_i)
   for j=1:n_samps
       if(round(j/emg_fsamp,4,'significant')==round(start_i(i),4,'significant'))
          muscle_time(i,1)=j;   
       end
       if(round(j/emg_fsamp,4,'significant')==round(end_i(i),4,'significant'))
          muscle_time(i,2)=j;
       end
   end
   muscle_time(i,3)=muscle_time(i,2)-muscle_time(i,1);
end

%% ricampionamento + andamento medio e variabilità

FILT_RES=struct;
samp=100;
SIGN_MED=struct;
for k=1:length(muscle_code)
    FILT_RES.(muscle_code{k})=[];
    SIGN_MED.(muscle_code{k})=[];
    res_muscle=[];
    res_muscle_m=[];
    if k<length(muscle_code)
        for i=1:length(muscle_time)
        tmp=[];
        res_tmp=[];
        tmp=FILT.(muscle_code{k})(1,muscle_time(i,1):muscle_time(i,2))';
        res_tmp=resample(tmp,samp,size(tmp,1))';
        res_muscle=[res_muscle res_tmp];
        res_muscle_m=[res_muscle_m;res_tmp]; % ogni riga rappresenta un ciclo passo, così le posso mediare tra loro
        end
        SIGN_MED.(muscle_code{k})(1,:)=mean(res_muscle_m);
        SIGN_MED.(muscle_code{k})(2,:)=std(res_muscle_m)/sqrt(length(SIGN_MED.(muscle_code{k})(1,:))); % errore standard
        FILT_RES.(muscle_code{k})(:,:)=res_muscle;
    else
        for i_ch=1:n_chs
            res_muscle=[];
            res_muscle_m=[];
            for i=1:length(muscle_time)
                tmp=[];
                res_tmp=[];
                tmp=FILT.(muscle_code{k})(i_ch,muscle_time(i,1):muscle_time(i,2))';
                res_tmp=resample(tmp,samp,size(tmp,1))';
                res_muscle=[res_muscle res_tmp];
                res_muscle_m=[res_muscle_m;res_tmp];
            end
            SIGN_MED.(muscle_code{k})(i_ch,:)=mean(res_muscle_m);
            SIGN_MED.(muscle_code{k})(i_ch+n_chs,:)=std(res_muscle_m)/sqrt(length(SIGN_MED.(muscle_code{k})(i_ch,:))); % errore standard
            FILT_RES.(muscle_code{k})(i_ch,:)=res_muscle;
        end
    end
end

%% plot segnali mediati
close all
t=1/emg_fsamp:1/emg_fsamp:length(SIGN_MED.BF)/1/emg_fsamp;
vnorm_in1=[];
vnorm_in=[];
vnorm_in=max(max(abs(SIGN_MED.GM)));
a=1;
figure()
for i=1:length(muscle_code)
    if i<length(muscle_code)
        subplot(3,3,i)
        SIGN_MED.(muscle_code{i})(1,:)=abs(SIGN_MED.(muscle_code{i})(1,:));
        plot(0:100/length(SIGN_MED.(muscle_code{i})(1,:)):100*(length(SIGN_MED.(muscle_code{i})(1,:))-1)/length(SIGN_MED.(muscle_code{i})(1,:)),(SIGN_MED.(muscle_code{i})(1,:)+SIGN_MED.(muscle_code{i})(2,:))*10^6,'b','LineWidth',1)
        hold on
        plot(0:100/length(SIGN_MED.(muscle_code{i})(1,:)):100*(length(SIGN_MED.(muscle_code{i})(1,:))-1)/length(SIGN_MED.(muscle_code{i})(1,:)),(SIGN_MED.(muscle_code{i})(1,:)-SIGN_MED.(muscle_code{i})(2,:))*10^6,'b','LineWidth',1)
        hold on   
        plot(0:100/length(SIGN_MED.(muscle_code{i})(1,:)):100*(length(SIGN_MED.(muscle_code{i})(1,:))-1)/length(SIGN_MED.(muscle_code{i})(1,:)),SIGN_MED.(muscle_code{i})(1,:)*10^6,'r','LineWidth',1.5), grid on,
        hold off
        title(['Segnale ', muscle_code{i}],' mediato'),xlabel('% of Gait'), ylabel('Attiv. muscolare (\muV)'),
        hold on

    else
       legend("errore standard sup","errore standard inf","segnale",'EdgeColor',[0 0 1]);
       for i_ch=1:n_chs
            
            if(i_ch==1||i_ch==2||i_ch==n_chs-1||i_ch==n_chs)
               subplot(3,3,i+a-1)
               a=a+1;
               SIGN_MED.(muscle_code{i})(i_ch,:)=abs(SIGN_MED.(muscle_code{i})(i_ch,:)); 
               plot(0:100/length(SIGN_MED.(muscle_code{i})(i_ch,:)):100*(length(SIGN_MED.(muscle_code{i})(i_ch,:))-1)/length(SIGN_MED.(muscle_code{i})(i_ch,:)),(SIGN_MED.(muscle_code{i})(i_ch,:)+SIGN_MED.(muscle_code{i})(i_ch+n_chs,:))*10^6,'b','LineWidth',1)
               hold on
               plot(0:100/length(SIGN_MED.(muscle_code{i})(i_ch,:)):100*(length(SIGN_MED.(muscle_code{i})(i_ch,:))-1)/length(SIGN_MED.(muscle_code{i})(i_ch,:)),(SIGN_MED.(muscle_code{i})(i_ch,:)-SIGN_MED.(muscle_code{i})(i_ch+n_chs,:))*10^6,'b','LineWidth',1)
               hold on   
               plot(0:100/length(SIGN_MED.(muscle_code{i})(i_ch,:)):100*(length(SIGN_MED.(muscle_code{i})(i_ch,:))-1)/length(SIGN_MED.(muscle_code{i})(i_ch,:)),SIGN_MED.(muscle_code{i})(i_ch,:)*10^6,'r','LineWidth',1.5), grid on
               hold off
               if i_ch<=2
                   txt=[muscle_code{i},' prossimale mediato al canale ',num2str(i_ch)]; 
                   xlabel('% of Gait'), ylabel('Attiv. muscolare (\muV)')
                   title(txt)
               else
                   txt=[muscle_code{i},' distale mediato al canale ',num2str(i_ch)]; 
                   xlabel('% of Gait'), ylabel('Attiv. muscolare (\muV)')
                   title(txt)
               end
            end
        end
    end
end

% Plot GM distale e prossimale
figure()
subplot(2,1,1)
SIGN_MED.(muscle_code{i})(1,:)=abs(SIGN_MED.(muscle_code{6})(1,:)); 
plot(0:100/length(SIGN_MED.(muscle_code{6})(1,:)):100*(length(SIGN_MED.(muscle_code{6})(1,:))-1)/length(SIGN_MED.(muscle_code{6})(1,:)),(SIGN_MED.(muscle_code{6})(1,:)+SIGN_MED.(muscle_code{6})(1+15,:))*10^6,'b','LineWidth',1)
hold on
plot(0:100/length(SIGN_MED.(muscle_code{6})(1,:)):100*(length(SIGN_MED.(muscle_code{6})(1,:))-1)/length(SIGN_MED.(muscle_code{6})(1,:)),(SIGN_MED.(muscle_code{6})(1,:)-SIGN_MED.(muscle_code{6})(1+15,:))*10^6,'b','LineWidth',1)
hold on   
plot(0:100/length(SIGN_MED.(muscle_code{6})(1,:)):100*(length(SIGN_MED.(muscle_code{6})(1,:))-1)/length(SIGN_MED.(muscle_code{6})(1,:)),SIGN_MED.(muscle_code{6})(1,:)*10^6,'r','LineWidth',1.5), grid on
txt=[muscle_code{6},' prossimale mediato al canale 1']; 
xlabel('% of Gait'), ylabel('Attiv. muscolare (\muV)')
title(txt)
hold off

subplot(2,1,2)
SIGN_MED.(muscle_code{i})(1,:)=abs(SIGN_MED.(muscle_code{6})(1,:)); 
plot(0:100/length(SIGN_MED.(muscle_code{6})(15,:)):100*(length(SIGN_MED.(muscle_code{6})(15,:))-1)/length(SIGN_MED.(muscle_code{6})(15,:)),(SIGN_MED.(muscle_code{6})(15,:)+SIGN_MED.(muscle_code{6})(15+15,:))*10^6,'b','LineWidth',1)
hold on
plot(0:100/length(SIGN_MED.(muscle_code{6})(15,:)):100*(length(SIGN_MED.(muscle_code{6})(15,:))-1)/length(SIGN_MED.(muscle_code{6})(15,:)),(SIGN_MED.(muscle_code{6})(15,:)-SIGN_MED.(muscle_code{6})(15+15,:))*10^6,'b','LineWidth',1)
hold on   
plot(0:100/length(SIGN_MED.(muscle_code{6})(15,:)):100*(length(SIGN_MED.(muscle_code{6})(15,:))-1)/length(SIGN_MED.(muscle_code{6})(15,:)),SIGN_MED.(muscle_code{6})(15,:)*10^6,'r','LineWidth',1.5), grid on
txt=[muscle_code{6},' distale mediato al canale 15']; 
xlabel('% of Gait'), ylabel('Attiv. muscolare (\muV)')
title(txt)

hold off

