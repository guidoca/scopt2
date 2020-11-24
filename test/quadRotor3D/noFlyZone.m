function [Out] = noFlyZone(X,Pobs)
%NOFLYZONE Summary of this function goes here
%   Detailed explanation goes here

Out = vecnorm((X(1:3,:)-Pobs),2,1);

end

