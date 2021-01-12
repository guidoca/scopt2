classdef aerodynamicModel
    %aerodynamicModel class with information on vehicle aerodynamics 
    %
    %   Receives coefficent of drag (Cd) and lift (Cl) data for Mach number
    %   (M), and provides different methods for computing the coefficients
    %   at any Mach number:
    %       - none: assumes there is no drag/lift
    %
    %       - constant: assumes drag/lift is constant
    %
    %       - interpolate: linear interpolation of each data point
    %
    %       - subsonic-fit: assumes three data points for subsonic
    %           regime, and approximates with piecewise function which is
    %           constant for the first two points and a second order
    %           polynomial for the last two, so that the derivative is
    %           zero at the second point.
    
    properties
        machData
        cdData
        clData
        dragMethod % 'none', 'constant', 'interpolate' or 'subsonic-fit'
        liftMethod % 'none', 'constant', 'interpolate' or 'subsonic-fit'
        area
    end
    
    properties (Dependent)
        includeDrag
        includeLift
    end
        
    properties (Access = private)
        dragInterpolant
        dragFit
        liftInterpolant
        liftFit
    end
    
    methods
        function obj = aerodynamicModel(machData,cdData,dragMethod,clData,liftMethod,area)
            if nargin >= 1
                obj.machData = machData;
            else
                obj.machData = [0,0.2,1];
            end
            
            if nargin >= 2
                obj.cdData = cdData;
            else
                obj.cdData = [1.0,1.0,0.6];
            end
            
            if nargin >= 3
                obj.dragMethod = dragMethod;
            else
                obj.dragMethod = 'none';
            end
            
            if nargin >= 4
                obj.clData = clData;
            else
                obj.clData = [0.0,0.0,0.0];
            end
            
            if nargin >= 5
                obj.liftMethod = liftMethod;
            else
                obj.liftMethod = 'none';
            end
            
            if nargin >= 6
                obj.area = area;
            else
                obj.area = 10;
            end
        end
        
        function [Cd,dCd] = computeCd(obj,mach,angleOfAttack)
            switch obj.dragMethod
                case 'none'
                    Cd = 0;
                    dCd = 0;
                    
                case 'constant'
                    Cd = obj.cdData(1);
                    dCd = 0;
                    
                case 'interpolate'
                    Cd = obj.dragInterpolant(mach);
                    dCd = NaN; % ??
                    
                case 'subsonic-fit'
                    if mach < obj.machData(2)
                        Cd = obj.cdData(2);
                        dCd = 0;
                    else
						Cd  = obj.dragFit(1)*mach.^2+obj.dragFit(2)*mach+obj.dragFit(3);
                        dCd = 2* obj.dragFit(1)*mach+obj.dragFit(2);
                    end
            end
        end
        
        function [Cl,dCl] = computeCl(obj,mach,angleOfAttack)
            switch obj.dragMethod
                case 'none'
                    Cl = 0;
                    dCl = 0;
                    
                case 'constant'
                    Cl = obj.clData(1);
                    dCl = 0;
                    
                case 'interpolate'
                    Cl = obj.liftInterpolant(mach);
                    dCl = NaN; % ??
                    
                case 'subsonic-fit'
                    if mach < obj.machData(2)
                        Cl = obj.clData(2);
                        dCl = 0;
                    else
                        Cl  = obj.liftFit(1)*mach.^2+obj.liftFit(2)*mach+obj.liftFit(3);
                        dCl = 2* obj.liftFit(1)*mach+obj.liftFit(2);
                    end
            end
        end
        
        function obj = set.dragMethod(obj,dragMethod)
            switch dragMethod
                case 'none'
                    
                case 'constant'
                    
                case 'interpolate'
                    obj.dragInterpolant = griddedInterpolant(obj.machData,obj.cdData);
                    
                case 'subsonic-fit'
                    obj.dragFit(1) = (obj.cdData(3)-obj.cdData(2))/(obj.machData(3)-obj.machData(2))^2;
                    obj.dragFit(2) = -obj.machData(2)*2*obj.dragFit(1);
                    obj.dragFit(3) = obj.cdData(2)+obj.machData(2)^2*obj.dragFit(1);
                otherwise
                    error('Invalid aerodynamic drag coefficient method.');
            end
            obj.dragMethod = dragMethod;
        end
        
        function obj = set.liftMethod(obj,liftMethod)
            switch liftMethod
                case 'none'
                    
                case 'constant'
                    
                case 'interpolate'
                    obj.dragInterpolant = griddedInterpolant(obj.machData,obj.clData);
                    
                case 'subsonic-fit'
                    obj.liftFit(1) = (obj.clData(3)-obj.clData(2))/(obj.machData(3)-obj.machData(2))^2;
                    obj.liftFit(2) = -obj.machData(2)*2*obj.liftFit(1);
                    obj.liftFit(3) = obj.clData(2)+obj.machData(2)^2*obj.liftFit(1);
                otherwise
                    error('Invalid aerodynamic lift coefficient method.');
            end
            obj.liftMethod = liftMethod;
        end
        
        function includeDrag = get.includeDrag(obj)
            if strcmp(obj.dragMethod,'none')
                includeDrag = false;
            else
                includeDrag = true;
            end
        end
        
        function includeLift = get.includeLift(obj)
            if strcmp(obj.liftMethod,'none')
                includeLift = false;
            else
                includeLift = true;
            end
        end
    end
end

