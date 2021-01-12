function B = controlMatrixRocketLandingFull(t,X,U,BodyMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nodes    = size(X,2);

%% Retrieve Environment
% Earth Shape Model 
% Thrust Model 
cex       = BodyMap.ThrustModel.ceffVac;

%% Compute State Derivatives
B = zeros(size(X,1),size(U,1),Nodes);

switch BodyMap.DynamicModel.controlVariable
    case 'thrust'
        M        = X(7,:);
        for ii = 1:Nodes
            B(:,:,ii) = [...
                zeros(3,4);...
                eye(3)./repmat(M(ii),3,1),zeros(3,1);...
                zeros(1,3),-1/cex;...
                ];
        end
    case 'acceleration'
        for ii = 1:Nodes
            B(:,:,ii) = [...
                zeros(3,4);...
                eye(3),zeros(3,1);...
                zeros(1,3),-1/cex;...
                ];
        end
end


end

