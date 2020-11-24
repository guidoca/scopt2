function [row,col,val,nEntries] = controlMatrixSpacecraft(t,X,U,BodyMap,ii)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
planet  = BodyMap.CentralBody      ;
vehicle = BodyMap.ControlledBody  ;
nodes   = size(X,2);

%% Retrieve Environment

% Mass Model
Mass     = BodyMap.(vehicle).MassModel.mass ;
%%
row = repmat((4:6)',1,nodes);
col = repmat((1:3)',1,nodes);

nEntries= size(row,1);

val = zeros(nEntries,nodes);

val( 1:3,:) = 1/Mass*ones(3,nodes);
end

