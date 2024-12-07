function [loc_ref_r,loc_ref_l,RTOE, LTOE]= calc_references_foot(traj,JOINTS)
% FOOT
% Input:        - marker 3D coordinates 
%               - RTOE
%               - LTOE
%               - RAJC
%               - RKJC
%               - LAJC
%               - LKJC
% Output:       - right system of reference coordinates
%               - right toe center coordinates 
%               - left system of reference coordinates
%               - left toe center coordinates 


%% RIGHT
RTOE=traj.RTOE;  % right center 
% z direction
dir_z=JOINTS.RA-RTOE;
ez_r=makeunit(dir_z);
% y direction 
temp=JOINTS.RK-RTOE;
dir_y=cross(ez_r,temp);
ey_r=makeunit(dir_y);
% x direction
ex_r=cross(ey_r,ez_r);      

loc_ref_r= [ex_r ey_r ez_r];

%% LEFT
LTOE=traj.LTOE;  % left center 

% z direction
dir_z=JOINTS.LA-LTOE;
ez_l=makeunit(dir_z);
% y direction
temp=JOINTS.LK-JOINTS.LA;
dir_y=cross(ez_l,temp);
ey_l=makeunit(dir_y);
% x direction
ex_l=cross(ey_l,ez_l); 

loc_ref_l= [ex_l ey_l ez_l];


end
