function [A,A1] = stateMatrixSolutionRocketLandingFull(t,X,U,BodyMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nodes    = size(X,2);
dt       = diff(t);

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

% Atmosphere Model
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
    cex          = BodyMap.ThrustModel.ceffVac;
    mdotbp       = -pamb.*BodyMap.ThrustModel.Ae./cex;
    dpdh         = BodyMap.AtmosphereModel.get_dpdh(h);
    dmdotbpdh    = -dpdh.*BodyMap.ThrustModel.Ae./cex;
else
    mdotbp       = zeros(1,Nodes);
    dmdotbpdh    = zeros(1,Nodes);
end

% Aerodynamic
if BodyMap.AerodynamicModel.includeDrag
    Vmag     = sqrt(V(1,:).^2+V(2,:).^2+V(3,:).^2);
    IndexZero = find(Vmag~=0);
    invVmag = zeros(1,Nodes);
    invVmag(IndexZero)   = 1./Vmag(IndexZero);
    
    mach    = Vmag./soundSpeed;
    [cD,dcDdM]      = BodyMap.AerodynamicModel.computeCd(mach);
    Area    = BodyMap.AerodynamicModel.area;
    Drag    = - 0.5*Area*repmat(rho.*cD.*Vmag,3,1).*V;
    aD      = Drag.*repmat(Minv,3,1);
    dMda    = -Vmag./soundSpeed.^2;
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
A  = zeros(size(X,1),size(X,1),Nodes-1);
A1 = zeros(size(X,1),size(X,1),Nodes-1);

switch BodyMap.DynamicModel.controlVariable
    case 'thrust'
        dadM = -a.*repmat(Minv,3,1);
        if BodyMap.DynamicModel.exactOde==0
            for ii = 1:Nodes-1
                A(:,:,ii) = [...
                    dt(ii)^2/2*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)^2/2*centrifugal_term,dt(ii)*eye(3)+dt(ii)^2/2*centripetal_term + dt(ii)^2/2*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)^2/2*daDdM(:,ii);...
                    dt(ii)*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)*centrifugal_term,dt(ii)*centripetal_term + dt(ii)*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)*daDdM(:,ii);...
                    dt(ii)*[dmdotbpdh(ii),zeros(1,2),zeros(1,4)] ;...
                    ];
            end
        elseif BodyMap.DynamicModel.exactOde==1
            for ii = 1:Nodes-1
                A(:,:,ii) = [...
                    dt(ii)^2/3*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)^2/3*centrifugal_term,dt(ii)*eye(3)+dt(ii)^2/3*centripetal_term + dt(ii)^2/3*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)^2/3*daDdM(:,ii);...
                    dt(ii)/2*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)/2*centrifugal_term,dt(ii)/2*centripetal_term + dt(ii)/2*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)/2*daDdM(:,ii);...
                    dt(ii)/2*[dmdotbpdh(ii),zeros(1,2),zeros(1,4)] ;...
                    ];
                A1(:,:,ii) = [...
                    dt(ii)^2/6*[ daDdh(:,ii+1) , zeros(3,2) ] + dt(ii)^2/6*centrifugal_term,dt(ii)^2/6*centripetal_term + dt(ii)^2/6*[daDxdV(:,ii+1)';daDydV(:,ii+1)';daDzdV(:,ii+1)'],dt(ii)^2/6*daDdM(:,ii+1);...
                    dt(ii)/2*[ daDdh(:,ii+1) , zeros(3,2) ] + dt(ii)/2*centrifugal_term,dt(ii)/2*centripetal_term + dt(ii)/2*[daDxdV(:,ii+1)';daDydV(:,ii+1)';daDzdV(:,ii+1)'],dt(ii)/2*daDdM(:,ii+1);...
                    dt(ii)/2*[dmdotbpdh(ii+1),zeros(1,2),zeros(1,4)] ;...
                    ];
            end
        else
            error('Dynamics Error: exactOde type not defined \n')
        end
        
        switch BodyMap.ThrustModel.collocationMethod
            case {'zoh','floor'}
                for ii = 1:Nodes-1
                    A(1:3,7,ii) = A(1:3,7,ii) + dt(ii)^2*dadM(:,ii);
                    A(4:6,7,ii) = A(4:6,7,ii) + dt(ii)*dadM(:,ii);
                end
            case {'foh','linear'}
                for ii = 1:Nodes-1
                    A(1:3,7,ii) = A(1:3,7,ii) + dt(ii)^2*dadM(:,ii)/3;
                    A(4:6,7,ii) = A(4:6,7,ii) + dt(ii)*dadM(:,ii)/2;
                    A1(1:3,7,ii) = A1(1:3,7,ii) + dt(ii)^2*dadM(:,ii+1)/6;
                    A1(4:6,7,ii) = A1(4:6,7,ii) + dt(ii)*dadM(:,ii+1)/2;
                end
        end
        
    case 'acceleration'
        if BodyMap.DynamicModel.exactOde==0
            for ii = 1:Nodes-1
                A(:,:,ii) = [...
                    dt(ii)^2/2*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)^2/2*centrifugal_term,dt(ii)*eye(3)+dt(ii)^2/2*centripetal_term + dt(ii)^2/2*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)^2/2*daDdZ(:,ii);...
                    dt(ii)*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)*centrifugal_term,dt(ii)*centripetal_term + dt(ii)*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)*daDdZ(:,ii);...
                    dt(ii)*[dmdotbpdh(ii)*Minv(ii),zeros(1,2),zeros(1,3),-mdotbp(ii)*Minv(ii)] ;...
                    ];
            end
        elseif BodyMap.DynamicModel.exactOde==1
            for ii = 1:Nodes-1
                A(:,:,ii) = [...
                    dt(ii)^2/3*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)^2/3*centrifugal_term,dt(ii)*eye(3)+dt(ii)^2/3*centripetal_term + dt(ii)^2/3*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)^2/3*daDdZ(:,ii);...
                    dt(ii)/2*[ daDdh(:,ii) , zeros(3,2) ] + dt(ii)/2*centrifugal_term,dt(ii)/2*centripetal_term + dt(ii)/2*[daDxdV(:,ii)';daDydV(:,ii)';daDzdV(:,ii)'],dt(ii)/2*daDdZ(:,ii);...
                    dt(ii)/2*[dmdotbpdh(ii)*Minv(ii),zeros(1,2),zeros(1,3),-mdotbp(ii)*Minv(ii)] ;...
                    ];
                A1(:,:,ii) = [...
                    dt(ii)^2/6*[ daDdh(:,ii+1) , zeros(3,2) ] + dt(ii)^2/6*centrifugal_term,dt(ii)^2/6*centripetal_term + dt(ii)^2/6*[daDxdV(:,ii+1)';daDydV(:,ii+1)';daDzdV(:,ii+1)'],dt(ii)^2/6*daDdZ(:,ii+1);...
                    dt(ii)/2*[ daDdh(:,ii+1) , zeros(3,2) ] + dt(ii)/2*centrifugal_term,dt(ii)/2*centripetal_term + dt(ii)/2*[daDxdV(:,ii+1)';daDydV(:,ii+1)';daDzdV(:,ii+1)'],dt(ii)/2*daDdZ(:,ii+1);...
                    dt(ii)/2*[dmdotbpdh(ii+1)*Minv(ii+1),zeros(1,2),zeros(1,3),-mdotbp(ii+1)*Minv(ii+1)] ;...
                    ];
            end
        else
            error('Dynamics Error: exactOde type not defined \n')
        end
        
end
end

