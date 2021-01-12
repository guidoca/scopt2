function A = stateMatrixRocketLandingFull(t,X,U,BodyMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nodes    = size(X,2);

R        = X(1:3,:);
V        = X(4:6,:);

switch BodyMap.DynamicModel.controlVariable
    case 'thrust' 
        M        = X(7,:);
        Minv     = 1./M;
        T        = U(1:3,:);
        Tmag     = U(4,:); 
        a        = T.*Minv; 
    case 'acceleration' 
        Z        = X(7,:);
        Minv     = exp(-Z);
        amag     = U(4,:); 
end

%% Retrieve Environment
% Earth Shape Model
h = R(1,:);

% Aerodynamic Model
if BodyMap.ThrustModel.includeBackPressure || BodyMap.AerodynamicModel.includeDrag
    [rho , pamb, ~ , soundSpeed] = BodyMap.AtmosphereModel.get_properties(h,[],[],[]);
end

% Orientation Model
if BodyMap.OrientationModel.includeCentrifugal
    centrifugal_term = - BodyMap.OrientationModel.rotRateSkew*BodyMap.OrientationModel.rotRateSkew;
else
    centrifugal_term = zeros(3,3) ;
end
if BodyMap.OrientationModel.includeCentripetal
    centripetal_term = - 2 * BodyMap.OrientationModel.rotRateSkew;
else
    centripetal_term = zeros(3,3) ;
end

% Thrust Model
if BodyMap.ThrustModel.includeBackPressure
    % get pressure derivatives with altitude
    cex         = BodyMap.ThrustModel.ceffVac; 
    dpdh         = BodyMap.AtmosphereModel.get_dpdh(h);
    mdotbp       = -pamb.*BodyMap.ThrustModel.Ae./cex;
    
    switch BodyMap.DynamicModel.controlVariable
        case 'thrust'
            dmdotbpdh    = -dpdh.*BodyMap.ThrustModel.Ae./cex;
        case 'acceleration'
            dmdotbpdh    = -dpdh.*BodyMap.ThrustModel.Ae./cex.*Minv;
    end
else
	dmdotbpdh    = zeros(1,Nodes);
	mdotbp    = zeros(1,Nodes);
end

% Aerodynamic
if BodyMap.AerodynamicModel.includeDrag
    Vmag     = sqrt(V(1,:).^2+V(2,:).^2+V(3,:).^2);
    IndexZero = find(Vmag~=0);
    invVmag = zeros(1,Nodes);
    invVmag(IndexZero)   = 1./Vmag(IndexZero);
    
    mach    = Vmag./soundSpeed;
    dMda    = -Vmag./soundSpeed.^2;
    [cD,dcDdM]      = BodyMap.AerodynamicModel.computeCd(mach);
    Area    = BodyMap.AerodynamicModel.area;
    Drag    = - 0.5*Area*repmat(rho.*cD.*Vmag,3,1).*V;
    aD      = Drag.*repmat(Minv,3,1);
    drhodh  = BodyMap.AtmosphereModel.get_drhodh(h);
    dadh    = BodyMap.AtmosphereModel.get_dadh(h);
    dcDdh   = dcDdM.*dMda.*dadh;
    dDdh    = - 0.5*Area*repmat((drhodh.*cD+rho.*dcDdh).*Vmag,3,1).*V;
    dDxdV   = - 0.5*Area*( V.*repmat( rho.*cD.*V(1,:).*invVmag ,3,1) + [Vmag;zeros(2,Nodes)]) - 0.5*Area*repmat(rho.*Vmag.*V(1,:).*dcDdM,3,1).*V.*repmat(invVmag./soundSpeed ,3,1);
    dDydV   = - 0.5*Area*( V.*repmat( rho.*cD.*V(2,:).*invVmag ,3,1) + [zeros(1,Nodes);Vmag;zeros(1,Nodes)]) - 0.5*Area*repmat(rho.*Vmag.*V(2,:).*dcDdM,3,1).*V.*repmat(invVmag./soundSpeed ,3,1);
    dDzdV   = - 0.5*Area*( V.*repmat( rho.*cD.*V(3,:).*invVmag ,3,1) + [zeros(2,Nodes);Vmag]) - 0.5*Area*repmat(rho.*Vmag.*V(3,:).*dcDdM,3,1).*V.*repmat(invVmag./soundSpeed ,3,1);

    switch BodyMap.DynamicModel.controlVariable
        case 'thrust'
            daDdM   = -aD.*repmat(Minv,3,1);
        case 'acceleration'
            daDdZ   = -aD;
    end
    daDdh   = dDdh .*Minv;
    daDxdV  = dDxdV.*Minv;
    daDydV  = dDydV.*Minv;
    daDzdV  = dDzdV.*Minv;
else 
    daDdM   = zeros(3,Nodes);
    daDdZ   = zeros(3,Nodes);
    daDdh   = zeros(3,Nodes);
    daDxdV  = zeros(3,Nodes);
    daDydV  = zeros(3,Nodes);
    daDzdV  = zeros(3,Nodes);
end
 
%% Compute State Derivatives
A = zeros(size(X,1),size(X,1),Nodes);
switch BodyMap.DynamicModel.controlVariable
    case 'thrust'
        dadM = -a.*repmat(Minv,3,1);
        for ii = 1:Nodes
            A(:,:,ii) = [...
                zeros(3),eye(3),zeros(3,1);...
                [daDdh(:,ii),zeros(3,2) ] + centrifugal_term,centripetal_term + [daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],daDdM(:,ii)+dadM(:,ii);...
                [dmdotbpdh(ii),zeros(1,2),zeros(1,4)] ;...
                ];
        end
    case 'acceleration'
        for ii = 1:Nodes
            A(:,:,ii) = [...
                zeros(3),eye(3),zeros(3,1);...
                [daDdh(:,ii),zeros(3,2) ] + centrifugal_term,centripetal_term + [daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],daDdZ(:,ii);...
                [dmdotbpdh(ii),zeros(1,2),zeros(1,3),-mdotbp(ii)*Minv(ii)] ;...
                ];
        end
end


end

