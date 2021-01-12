classdef orientationModel
    % Orientation model for coordinate system choice
    
    properties (Constant)
        rotRate_E = 7.2921159e-5 ; % rad/s
    end
    properties
        rotRate = 7.2921159e-5 ;
        rotRateDirection = [ 0 ; 1 ; 0] ; 
		
		includeCentrifugal = 0 ;
		includeCentripetal = 0 ;
    end
    properties (Dependent)
        rotRateVector 
        rotRateSkew
    end
    
    methods 
		function [a] = computeRotatingAcceleration(obj,R,V)  
			a   = zeros(size(R)) ;
			if obj.includeCentrifugal
				a = a + computeCentrifugalAcceleration(obj,R);
			end
			if obj.includeCentripetal
				a = a + computeCentripetalAcceleration(obj,V);
			end
		end	
		function [a] = computeCentrifugalAcceleration(obj,R) 
			a = - obj.rotRateSkew*obj.rotRateSkew*R;
		end		 
		function [a] = computeCentripetalAcceleration(obj,V) 
            a = - 2 * obj.rotRateSkew * V;
		end		 
        function rotRateVector = get.rotRateVector(obj)
            rotRateVector = obj.rotRateDirection*obj.rotRate;
        end
        function rotRateSkew = get.rotRateSkew(obj)
            rotRateSkew = [...
                0 -obj.rotRateVector(3) obj.rotRateVector(2);...
                obj.rotRateVector(3) 0 -obj.rotRateVector(1);...
                -obj.rotRateVector(2) obj.rotRateVector(1) 0] ;
        end
    end
    
end

