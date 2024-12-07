function [loc_ref, JC,RHJC_g,LHJC_g]= calc_references_pelvis(traj) 

% PELVIS
% Input:        - marker 3D coordinates 
% Output:       - system of reference coordinates
%               - center coordinates 
%               - RHJC (right hip joint centre coordinates)
%               - LHJC (left hip joint centre coordinates)



%PELVIS
origin=(traj.RASI+traj.LASI)/2;
% y diretion 
dir_y=traj.LASI- origin; 
ey=makeunit(dir_y);             

% z diretion
sacro=((traj.RPSI)+(traj.LPSI))./2;
temp=sacro-origin;
dir_z=cross(ey,temp);
ez=makeunit(dir_z);

% x diretion
ex=cross(ey,ez);
loc_ref= [ex ey ez];
 
% hip joint center left (misure in mm)
%Antropometric masures (in mm) for hip joint center estimation;
Antro.LASI_RASI_dist= 240; %distance between LASI and RASI
Antro.leg_length= 980; %leg length

beta=0.314;
teta=0.5;
atd=0.1288*980-48.56;
c=Antro.leg_length*0.115-15.3;
aa=Antro.LASI_RASI_dist/2;
mm= 10/2; 

LHJCx=c*cos(teta)*sin(beta)-(atd+mm)*cos(beta);
LHJCy=-(c*sin(teta)-aa);
LHJCz=-c*cos(teta)*cos(beta)-(atd+mm)*sin(beta);

LHJC=[LHJCx LHJCy LHJCz];
LHJC_g=LHJC+origin;

% for the right joint centre, the y  is negated
RHJCy=(c*sin(teta)-aa);              
RHJC=[LHJCx RHJCy LHJCz];
RHJC_g=RHJC+origin;

JC=(LHJC_g+RHJC_g)./2;

end
