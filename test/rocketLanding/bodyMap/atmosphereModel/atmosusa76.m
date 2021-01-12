function [T, a, P, Rho,derivatives] = atmosusa76( height ,flagDerivative)
%atmosusa76 Retrieve atmospheric properties based on US76 standard
%atmosphere. Currently only working up to 11 km
%   Detailed explanation goes here 
if ~exist('flagDerivative','var')
    flagDerivative = 0;
end
derivatives = struct;
%% Constants and Parameters
R0e   = 6356.766e3; % m Equatorial Altitude for US76 Atmosphere Model
g0    = 9.80665  ; % m/s^2 Sea level Surface Gravity
M0    = computeMolarMass(0)  ; % %kg/kmol Temporary, Depends on Altitude for unmixed layers
Rstar = 8.31432  ;% J/mol/K Universal Gas Constant 
gamma = 1.4      ; % Ratio of specific heat of air at constant pressure
T_c   = 263.1905 ; % K
A     = -76.3232 ;
b     = 19.9428  ;

% Table
I_i   = [0      ]  ; % Layer Index
H_i   = [0      ]  ; % m Altitude
TM_i  = [288.15 ]  ; % K Layer Base Temperature
L_i   = [-6.5e-3]  ; % K/m % Layer Lapse Rate
P_i   = [1.013e5]  ; % Pa % Layer Base Pressure
%%
% Compute Geopotential Altitude
H   = R0e.*height./(R0e+height);
if flagDerivative % Derivative
    dHdh = R0e./(R0e+height) - R0e.*height./(R0e+height).^2 ;  
end
% Compute Molecular Scale Temperature
TM  = TM_i(1) + L_i(1)*(H-H_i(1));
if flagDerivative  % Derivative
    dTMdh = L_i(1).*dHdh;
end
% Compute Molar Mass
M   = computeMolarMass(H);
if flagDerivative % Derivative
    dMdH = computeMolarMassDerivative(H) ; 
end
% Compute Temperature
T   = TM .* M ./ M0 ;
if flagDerivative % Derivative
    dTdh = ( dTMdh .* M + TM .* dMdH .* dHdh ) ./ M0 ;
    derivatives.dTdh = dTdh;
end
% Compute Pressure
P   = P_i.*exp(-g0.*M0.*(H - H_i(1))/Rstar./TM) ;
if flagDerivative % Derivative
    dPdh = -g0.*M0/Rstar./TM.*P_i.*exp(-g0.*M0.*(H - H_i(1))/Rstar./TM).*(dHdh - (H - H_i(1))./TM.^2.*dTMdh) ;
    derivatives.dPdh = dPdh;
end
% Compute Density
Rho = M0 / Rstar * P ./ TM  ;
if flagDerivative % Derivative
    dRhodh = M0 / Rstar *(dPdh./TM - P ./TM.^2.*dTMdh);
    derivatives.dRhodh = dRhodh;
end
% Compute Speed of Sound
a   = sqrt(gamma*Rstar/M0.*TM);
if flagDerivative % Derivative
    dadh = sqrt(gamma*Rstar/M0./TM).*1/2.*dTMdh ; 
    derivatives.dadh = dadh;
end
end

function M = computeMolarMass(H)
% Computes molar mass based on altitude and tabluated fraction table
% currently only retrieves sea level one
M0  = 28.9644/1000; % %kg/kmol  
    
M   = M0*ones(size(H));
end

function dMdH = computeMolarMassDerivative(H)
% Computes molar mass derivative based on altitude and tabluated fraction table
% currently only retrieves sea level one
 
dMdH   = zeros(size(H));
end