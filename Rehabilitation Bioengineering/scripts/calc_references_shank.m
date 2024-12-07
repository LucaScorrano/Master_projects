function [loc_ref_r,loc_ref_rnt,loc_ref_l,loc_ref_lnt,RA,LA]= calc_references_shank(traj,JOINTS,femur_r,femur_l)
% THIGH
% Input:        - marker 3D coordinates 
%               - RKJC
%               - LKJC
%               - right system of reference coordinates
%               - leftt system of reference coordinates
% Output:       - right torisoned ankle system of reference coordinates
%               - right untorisoned ankle system of reference coordinates
%               - left torisoned ankle system of reference coordinates
%               - left untorisoned ankle system of reference coordinates
%               - right ankle centre coordinates 
%               - left ankle centre coordinates 


%Antropoometric masures (in mm) for ankle joint center estimation
Antro.ankle_width= 80; %ankle width in mm
mm=10/2;
AOS=mm+Antro.ankle_width/2;
%% Right
a=traj.RTIB;  
b=JOINTS.RK;    
c=traj.RANK;      
RA=chordPiG (a,b,c,AOS);
% torsionata
% z direction
dir_z=JOINTS.RK-RA;
ez_r=makeunit(dir_z);
% x direction
temp=traj.RTIB-RA;
dir_x=cross(ez_r,temp);
ex_r=makeunit(dir_x);
% y direction
ey_r=cross(ez_r,ex_r);

loc_ref_r= [ex_r ey_r ez_r];

% non torsionata

% z direction
dir_z=JOINTS.RK-RA;
e_zr=makeunit(dir_z);
% x direction
e_xr=femur_r(:,1:3);
% y direction
e_yr=cross(e_zr,e_xr);

loc_ref_rnt= [e_xr e_yr e_zr];

%% Left
a=traj.LTIB;  
b=JOINTS.LK;    
c=traj.LANK;  

LA=chordPiG (a,b,c,AOS);

% torsionata

% z direction
dir_z=JOINTS.LK-LA;
ez_l=makeunit(dir_z);
% x direction
%temp=LA-traj.LANK;
temp=traj.LTIB-LA;
dir_x=cross(temp,ez_l);
ex_l=makeunit(dir_x);
% y direction
ey_l=cross(ez_l,ex_l);

loc_ref_lnt= [ex_l  ey_l ez_l];

% non torsionata

% z direction
dir_z=JOINTS.LK-LA;
e_zl=makeunit(dir_z);
% x direction
e_xl=femur_l(:,1:3);
% y direction
e_yl=cross(e_zl,e_xl);

loc_ref_l= [e_xl e_yl e_zl];

end
 