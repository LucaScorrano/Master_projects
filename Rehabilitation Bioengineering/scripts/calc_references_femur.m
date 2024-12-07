 function [loc_ref_r, loc_ref_l,RK, LK]= calc_references_femur(traj,JOINTS)

% FEMUR
% Input:        - marker 3D coordinates 
%               - RHJC
%               - LHJC
% Output:       - right system of reference coordinates
%               - leftt system of reference coordinates
%               - right knee center coordinates 
%               - left knee center coordinates 


% RIGHT
%Antropoometric masures (in mm) for knee joint center estimation
Antro.knee_width= 100; %knee width in mm
mm=10/2;
delta=Antro.knee_width/2+mm;
a=traj.RTHI;         % wand marker 
b=JOINTS.RH;         % marker distale
c=traj.RKNE;         % marker prossimale
RK=chordPiG (a,b,c,delta);     

% z direction
dir_z=JOINTS.RH-RK;
ez_r=makeunit(dir_z);
% x direction
temp=traj.RTHI-traj.RKNE;
dir_x=cross(ez_r,temp);
ex_r=makeunit(dir_x);
% y direction
ey_r=cross(ez_r,ex_r);

loc_ref_r= [ex_r ey_r ez_r];

%% LEFT
a=traj.LTHI;         % wand marker 
b=JOINTS.LH;         % marker distale
c=traj.LKNE;         % marker prossimale
delta=100/2+10/2;
LK=chordPiG (a,b,c,delta);     


% z direction
dir_z=JOINTS.LH-LK;
ez_l=makeunit(dir_z);
% x direction
temp=traj.LTHI-traj.LKNE;
dir_x=cross(temp,ez_l);
ex_l=makeunit(dir_x);
% y direction
ey_l=cross(ez_l,ex_l);

loc_ref_l= [ex_l ey_l ez_l];
 end

