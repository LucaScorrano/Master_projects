function [angles]= joint_angles(proxi,dist)
%   INPUT:  - matrice di rotazione prossimale
%           - matrice di rotazione distale
%   OUTPUT: - angoli articolari


   R=[];
   R=proxi'*dist;
   alpha=asind(R(3,2));
   beta=asind((-R(3,1))/cosd(alpha));
   gamma=asind((-R(1,2))/cosd(alpha));

   angles=[alpha' beta' gamma'];

end

