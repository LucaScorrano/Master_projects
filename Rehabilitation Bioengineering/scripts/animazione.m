clear all 
close all 
clc 
 
%%define global variables to reduce parameters in the functions 
global markers  
global joints 
global sides 
global ref
%% Loading Data 
fsamp=100; %sample frequency 
markers= {'RASI'; 'LASI'; 'RPSI'; 'LPSI';'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'}; 
joints={'PC'; 'RK'; 'LK'; 'RA'; 'LA';'RTOE';'LTOE'}; 
sides= {'left'; 'right'}; 
ref={'pelvis','femur_r','femur_l','thigh_r','thigh_l','foot_r','foot_l'};

%Load marker coordinates 
[traj]=load_kin_data(fsamp); 
%swap columns 1 and 2 in order to have medio lateral axis along y 
for i_m=1:length(markers) 
    traj.(markers{i_m})= traj.(markers{i_m})(:,[2 1 3]); 
end 
 
%% Animazione Iniziale 
 
%suddivisione in x,y,z per ciascun marker 
for i=1:length(markers) 
x(:,i)=traj.(markers{i})(:,1); 
y(:,i)=traj.(markers{i})(:,2); 
z(:,i)=traj.(markers{i})(:,3); 
end  

%definizione vettori per unione marker con linee 
lin1=[1 1 2 3 1 5 6 7 8 9 10 2 11 12 13 14 15 16]; 
lin2=[2 3 4 4 5 6 7 8 9 10 8 11 12 13 14 15 16 14]; 
 
JOINTS=struct;
loc_ref=struct;
[loc_ref.pelvis,JOINTS.PC,JOINTS.RH,JOINTS.LH]=calc_references_pelvis(traj);
[loc_ref.femur_r,loc_ref.femur_l,JOINTS.RK,JOINTS.LK]=calc_references_femur(traj,JOINTS);
[loc_ref.thigh_r, loc_ref.thigh_r_nt,loc_ref.thigh_l, loc_ref.thigh_l_nt,JOINTS.RA, JOINTS.LA]= calc_references_shank(traj,JOINTS,loc_ref.femur_r,loc_ref.femur_l);
[loc_ref.foot_r,loc_ref.foot_l,JOINTS.RTOE ,JOINTS.LTOE]= calc_references_foot(traj,JOINTS);

JOINTS= rmfield(JOINTS,'RH');
JOINTS=rmfield(JOINTS,'LH');
loc_ref=rmfield(loc_ref,'thigh_r_nt');
loc_ref=rmfield(loc_ref,'thigh_l_nt');
%% Animazione con Sistemi di Riferimento 
 
% vettori per unire i marker 
lin1=[1 1 2 3 1 5 6 7 8 9 10 2 11 12 13 14 15 16]; 
lin2=[2 3 4 4 5 6 7 8 9 10 8 11 12 13 14 15 16 14]; 
for i=1:length(markers) 
    h_marker(i)=plot3(x(254,i),y(254,i),z(254,i),'ob'); 
    hold on     
    grid on  
end 
for i=1:length(lin1) 
    h_lin(i)=line([x(254,lin1(i)) x(254,lin2(i))],[y(254,lin1(i)) y(254,lin2(i))],[z(254,lin1(i)) z(254,lin2(i))],'LineStyle','--'); 
    hold on 
end 
hold on 
axis([-4000 6000 -2000 2000 -500 1500]) 
% asse x in rosso
% asse y in blu
% asse z in verde
for i=1:length(joints) 
    h1(i)=plot3(JOINTS.(joints{i})(1,1),JOINTS.(joints{i})(1,2),JOINTS.(joints{i})(1,3),'ok'); 
    hold on     
    h2(i)=quiver3(JOINTS.(joints{i})(1,1),JOINTS.(joints{i})(1,2),JOINTS.(joints{i})(1,3),loc_ref.(ref{i})(1,1),loc_ref.(ref{i})(1,4),loc_ref.(ref{i})(1,7),150,'r'); 
    hold on 
    h3(i)=quiver3(JOINTS.(joints{i})(1,1),JOINTS.(joints{i})(1,2),JOINTS.(joints{i})(1,3),loc_ref.(ref{i})(1,2),loc_ref.(ref{i})(1,5),loc_ref.(ref{i})(1,8),150,'b'); 
    hold on 
    h4(i)=quiver3(JOINTS.(joints{i})(1,1),JOINTS.(joints{i})(1,2),JOINTS.(joints{i})(1,3),loc_ref.(ref{i})(1,3),loc_ref.(ref{i})(1,6),loc_ref.(ref{i})(1,9),150,'g'); 
    hold on 
    grid on 
    legend("","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","centro","asse x","asse y", "asse z")
end 
hold on 
title('Sistemi di Riferimento'), xlabel('asse x (mm)'), ylabel('asse y (mm)'), zlabel('asse z (mm)')
axis([-4000 6000 -2000 2000 -500 1500]) 
for k=280:5:500 
for i=1:length(lin1) 
    if i<=length(markers) 
    set(h_marker(i),'XData',x(k,i),'YData',y(k,i),'ZData',z(k,i),'Marker','o','Color','b') 
    end 
    set(h_lin(i),'XData',[x(k,lin1(i)) x(k,lin2(i))],'YData',[y(k,lin1(i)) y(k,lin2(i))],'ZData',[z(k,lin1(i)) z(k,lin2(i))],'LineStyle','--') 
   
end 
for i=1:length(joints) 
    set(h1(i),'XData',JOINTS.(joints{i})(k,1),'YData',JOINTS.(joints{i})(k,2),'ZData',JOINTS.(joints{i})(k,3)) 
    set(h2(i),'XData',JOINTS.(joints{i})(k,1),'YData',JOINTS.(joints{i})(k,2),'ZData',JOINTS.(joints{i})(k,3),'UData',loc_ref.(ref{i})(k,1),'VData',loc_ref.(ref{i})(k,2),'WData',loc_ref.(ref{i})(k,3)); % asse x 
    set(h3(i),'XData',JOINTS.(joints{i})(k,1),'YData',JOINTS.(joints{i})(k,2),'ZData',JOINTS.(joints{i})(k,3),'UData',loc_ref.(ref{i})(k,4),'VData',loc_ref.(ref{i})(k,5),'WData',loc_ref.(ref{i})(k,6)); % asse y
    set(h4(i),'XData',JOINTS.(joints{i})(k,1),'YData',JOINTS.(joints{i})(k,2),'ZData',JOINTS.(joints{i})(k,3),'UData',loc_ref.(ref{i})(k,7),'VData',loc_ref.(ref{i})(k,8),'WData',loc_ref.(ref{i})(k,9)); % asse z
    
end 
pause(0.1)  
end





