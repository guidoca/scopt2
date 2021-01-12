classdef atmosphereModel 
    %ATMOSPHEREModel Model class with information atmospheric model properties
    %   Default is atmosphere with no altitude dependency
    %
    %   Receives the altitude, longitude and latitude and provides the atmospheric
    %   properties density, speed of sound, temperature and pressure:
    %       - const: constant atmosphere model with no dependency on altitude, longitude and latitude
    %
    %       - exp: Exponential atmospheric model with exponential dependency on geodetic altitude and constant temperature
    %
    %       - tabulated: Uses a tabulated atmosphere dependent on geodetic altitude. Default is the US76 Standaerd Atmosphere from ASTOS
    %
    %       - nrlmsise00: Uses matlab's nrlmsise00 reference atmosphere  with dependency on geodetic altitude, longitude and geodetic latitude
	
	% Can also provide the first and second derivatives with respect the dependent variables
     
    properties (Dependent)
        constantTemperature
        constantPressure
    end
    properties
        model = 'const';
        g0   = 9.807;
        rho0 = 1.225; % air density at sea level [kg/m^3]
        Hs   = 7.2*1e3;
        gasConstant = 287;
        gamma = 1.4;
        
        tabulatedAtmospherePath = '';
        interpmethod = 'spline';
        tabulatedAtmosphereModelPath = char;
        tempDatabase        = struct;
        database            = struct;
        database_derivative = struct;
        database_secondderivative = struct;
        
        constantTime = 1;
        year         = 2020;
        dayOfYear    = 1;
        UTseconds    = 0;
        
        densityInterpolant;
        speedOfSoundInterpolant;
        temperatureInterpolant;
        pressureInterpolant;
        
        densityError      = 0;
        speedOfSoundError = 0;
        temperatureError  = 0;
        pressureError     = 0;
        
        action = 'None'
        
        h_altitude = 1;
        
    end
    
    methods
        function obj = atmosphereModel(model,arg1) % Constructor
            if nargin == 0
                model = obj.model;
            end 
            switch model
                case 'const' % Independent Atmosphere
                case 'exp'
                case 'tabulated'
                    if nargin >= 2
                        obj.interpmethod =arg1;
                    end 
                case 'us76'
                case 'nrlmsise00'
                otherwise
                    error('Invalid Atmospheric Model');
            end
            
            obj.model    = model ;
        end
        
        function obj = set.model(obj,model)
            obj.model    = model ;
            switch model
                case 'const' % Independent Atmosphere
                case 'exp'
                case 'tabulated'
                    obj.tabulatedAtmosphereModelPath = fullfile(obj.tabulatedAtmospherePath, 'USSA1976Until130kmASTOS_HD.txt');
                    obj.database                     = load(obj.tabulatedAtmosphereModelPath);
                    obj = updateDatabaseInterpolantMethod(obj);
                case 'us76'
                    % From U.S. Standard Atmosphere, 1976: https://ntrs.nasa.gov/search.jsp?R=19770009539
                case 'nrlmsise00'
                otherwise
                    error('Invalid Atmospheric Model');
            end
        end
        function obj = set.interpmethod(obj,interpmethod)
            if strcmp('tabulated',obj.model)
                obj.interpmethod = interpmethod ; 
                obj = updateDatabaseInterpolantMethod(obj,interpmethod);
            else
                warning(['Interpolation method could not be set as current amtosphere model is not US76. Staying with ' obj.interpmethod] )
            end
        end
        
        
        function obj = updateDatabaseInterpolantMethod(obj,interpmethod)
            if nargin <2
                interpmethod = obj.interpmethod ; 
            end
            if ~strcmp(interpmethod,obj.interpmethod)
                obj.interpmethod = interpmethod ; 
            end
            obj.densityInterpolant      = griddedInterpolant(obj.database(:,1),obj.database(:,2),interpmethod);
            obj.pressureInterpolant     = griddedInterpolant(obj.database(:,1),obj.database(:,3),interpmethod);
            obj.temperatureInterpolant  = griddedInterpolant(obj.database(:,1),obj.database(:,4),interpmethod);
            obj.speedOfSoundInterpolant = griddedInterpolant(obj.database(:,1),obj.database(:,5),interpmethod);
        end
        function [rho , p, T , a] = get_properties(obj,h,lon,lat,time)
            % Altitude in meters 
            switch obj.model 
                case 'const' % Independent Atmosphere
                    rho = repmat(obj.rho0,size(h));
                    T   = repmat(obj.constantTemperature,size(h));
                    p   = rho.*T*obj.gasConstant;
                    a   = sqrt(obj.gamma*obj.gasConstant*T);
                case 'exp' % Exponential Atmosphere
                    h   = max(h,0);
                    rho = obj.rho0*exp(-h/obj.Hs);
                    T   = ones(size(h))*obj.constantTemperature;
                    p   = rho.*T*obj.gasConstant;
                    a   = sqrt(obj.gamma*obj.gasConstant*T);
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    h = min(max(h,min(obj.database(:,1))),max(obj.database(:,1))); 
                    rho= obj.densityInterpolant(h);
                    p  = obj.pressureInterpolant(h);
                    T  = obj.temperatureInterpolant(h);
                    a  = obj.speedOfSoundInterpolant(h);
                case 'us76'
                    [T, a, p, rho] = atmosusa76( h );
                case 'nrlmsise00' % nrlmsise00
                    if obj.constantTime
                        time.year      = obj.year*ones(size(h));
                        time.dayOfYear = obj.dayOfYear*ones(size(h));
                        time.UTseconds = obj.UTseconds*ones(size(h));
                    end
                    [TVec , rhoVec] = atmosnrlmsise00(h(:), rad2deg(lat(:)), rad2deg(lon(:)), time.year(:), time.dayOfYear(:), time.UTseconds(:),obj.action);
                    rho = reshape(rhoVec(:,6),size(h))                   ;
                    T   = reshape(TVec(:,2),size(h))                     ;
                    p   = rho.*T*obj.gasConstant            ;
                    a   = sqrt(obj.gamma*obj.gasConstant*T) ;
            end
            
            if obj.densityError ~= 0 
                rho = rho*(1+obj.densityError/100);
            end
            if obj.speedOfSoundError ~= 0 
                a = a*(1+obj.speedOfSoundError/100);
            end
        end
        
        function [drhodh , dpdh, dTdh , dadh] = get_derivatives(obj,h)
            % Altitude in meters 
            drhodh = get_drhodh(obj,h) ;
            dpdh   = get_dpdh(obj,h)   ;
            dTdh   = get_dTdh(obj,h)   ;
            dadh   = get_dadh(obj,h)   ;
             
        end
        
        function rho= get_rho(obj,h,lon,lat,time)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere
                    rho = repmat(obj.rho0,size(h)) ;
                case 'exp' % Exponential Atmosphere
                    rho = obj.rho0*exp(-max(h,0)/obj.Hs);
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    rho = obj.densityInterpolant(min(max(h,min(obj.database(:,1))),max(obj.database(:,1))));
                case 'us76'
                    [~, ~, ~, rho] = atmosusa76( h );
                case 'nrlmsise00'
                    if nargin == 5
                        rho = get_properties(obj,h,lon,lat,time);
                    elseif nargin == 4
                        rho = get_properties(obj,h,lon,lat);
                    end
            end
            if obj.densityError ~= 0 
                rho = rho*(1+obj.densityError/100);
            end
        end
        
        
        function speedOfSound = get_asound(obj,h,lon,lat,time)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere
                    speedOfSound = repmat(sqrt(obj.gamma*obj.gasConstant*obj.constantTemperature),size(h)) ;
                case 'exp' % Exponential Atmosphere
                    speedOfSound   = ones(size(h))*sqrt(obj.gamma*obj.gasConstant*obj.constantTemperature);
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    speedOfSound   = obj.speedOfSoundInterpolant(min(max(h,min(obj.database(:,1))),max(obj.database(:,1))));
                case 'us76'
                    [~, speedOfSound] = atmosusa76( h );
                case 'nrlmsise00'
                    if nargin == 5
                        [~ , ~, ~ , speedOfSound] = get_properties(obj,h,lon,lat,time);
                    elseif nargin == 4
                        [~ , ~, ~ , speedOfSound] = get_properties(obj,h,lon,lat);
                    end
            end
            if obj.speedOfSoundError ~= 0 
                speedOfSound = speedOfSound*(1+obj.speedOfSoundError/100);
            end
        end
        
        function p = get_p(obj,h,lon,lat,time)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere  
                    p   = repmat(obj.rho0.*obj.constantTemperature*obj.gasConstant,size(h));
                case 'exp' % Exponential Atmosphere
                    p = obj.rho0*exp(-max(h,0)/obj.Hs).*obj.constantTemperature*obj.gasConstant;
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    p = obj.pressureInterpolant(min(max(h,min(obj.database(:,1))),max(obj.database(:,1))));
                case 'us76'
                    [~, ~ , p] = atmosusa76( h );
                case 'nrlmsise00'
                    if nargin == 5
                        [~,p] = get_properties(obj,h,lon,lat,time);
                    elseif nargin == 4
                        [~,p] = get_properties(obj,h,lon,lat);
                    end
            end
            if obj.pressureError ~= 0 
                p = p*(1+obj.pressureError/100);
            end
        end
        
        
        function T = get_T(obj,h,lon,lat,time)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere
                    T = repmat(obj.constantTemperature,size(h)) ;
                case 'exp' % Exponential Atmosphere
                    T = ones(size(h))*obj.constantTemperature;
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    T = obj.temperatureInterpolant(min(max(h,min(obj.database(:,1))),max(obj.database(:,1))));
                case 'us76'
                    T = atmosusa76( h );
                case 'nrlmsise00'
                    if nargin == 5
                        [~ , ~, T] = get_properties(obj,h,lon,lat,time);
                    elseif nargin == 4
                        [~ , ~, T] = get_properties(obj,h,lon,lat);
                    end
            end
            if obj.temperatureError ~= 0 
                T = T*(1+obj.temperatureError/100);
            end
        end
        
        function machNumber = get_machNumber(obj,speed,h,lon,lat,time)
            if nargin <6
                machNumber   = speed./get_asound(obj,h,lon,lat);
            else
                machNumber   = speed./get_asound(obj,h,lon,lat,time);
            end
        end
        function speed = get_speed(obj,machNumber,h,lon,lat,time)
            if nargin <6
            speed   = machNumber.*get_asound(obj,h,lon,lat);
            else
            speed   = machNumber.*get_asound(obj,h,lon,lat,time);
            end
        end
        
        
        function drhodh = get_drhodh(obj,h)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere
					if exist('h','var')
						drhodh = zeros(size(h)) ;
					else
						drhodh = 0 ;
					end
                case 'exp' % Exponential Atmosphere
                    rho    = get_rho(obj,max(h,0));
                    drhodh = -1/obj.Hs*rho;
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    drhodh = (get_rho(obj,h + obj.h_altitude/2) - get_rho(obj,h - obj.h_altitude/2))/obj.h_altitude;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    drhodh = derivatives.dRhodh;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        
        function d2rhodh2 = get_d2rhodh2(obj,h)
            % Altitude in meters 
            switch obj.model  
                case 'const' % Independent Atmosphere
					if exist('h','var')
						d2rhodh2 = zeros(size(h)) ;
					else
						d2rhodh2 = 0;
					end
                case 'exp' % Exponential Atmosphere
                    rho      = get_rho(obj,max(h,0));
                    d2rhodh2 = 1/obj.Hs^2*rho;
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    d2rhodh2 = (get_rho(obj,h + obj.h_altitude) + get_rho(obj,h - obj.h_altitude) - 2*get_rho(obj,h))/obj.h_altitude^2;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    d2rhodh2 = derivatives.d2Rhodh2;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        function dTdh = get_dTdh(obj,h)
            % Altitude in meters
            switch obj.model
                case 'const' % Independent Atmosphere
                    if exist('h','var')
                        dTdh = zeros(size(h)) ;
                    else
                        dTdh = 0;
                    end
                case 'exp' % Exponential Atmosphere
                    dTdh = zeros(size(h));
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    dTdh = (get_T(obj,h + obj.h_altitude/2) - get_T(obj,h - obj.h_altitude/2))/obj.h_altitude;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    dTdh = derivatives.dTdh;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        function dadh = get_dadh(obj,h)
            % Altitude in meters
            switch obj.model
                case 'const' % Independent Atmosphere
                    if exist('h','var')
                        dadh = zeros(size(h)) ;
                    else
                        dadh = 0;
                    end
                case 'exp' % Exponential Atmosphere
                    dadh = zeros(size(h));
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    dadh = (get_asound(obj,h + obj.h_altitude/2) - get_asound(obj,h - obj.h_altitude/2))/obj.h_altitude;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    dadh = derivatives.dadh;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        function d2adh2 = get_d2adh2(obj,h)
            % Altitude in meters 
            switch obj.model
				case 'const' % Independent Atmosphere
					if exist('h','var')
						d2adh2 = zeros(size(h)) ;
					else
						d2adh2 = 0;
					end	
                case 'exp' % Exponential Atmosphere
                    d2adh2 =  zeros(size(h));
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    d2adh2 = (get_asound(obj,h + obj.h_altitude) + get_asound(obj,h - obj.h_altitude) - 2*get_asound(obj,h))/obj.h_altitude^2;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    d2adh2 = derivatives.d2adh2;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        function dpdh = get_dpdh(obj,h)
            % Altitude in meters  
            switch obj.model
                case 'const' % Independent Atmosphere
					if exist('h','var')
						dpdh = zeros(size(h)) ;
					else
						dpdh = 0 ;
					end
                case 'exp' % Exponential Atmosphere
                    dpdh = obj.get_drhodh(h)*obj.constantTemperature*obj.gasConstant;
                case 'tabulated' % Tabulated Atmosphere Model. Default is US76
                    dpdh = (get_p(obj,h + obj.h_altitude/2) - get_p(obj,h - obj.h_altitude/2))/obj.h_altitude;
                case 'us76'
                    [~,~,~,~,derivatives] = atmosusa76( h ,1);
                    dpdh = derivatives.dPdh;
                otherwise
                    error('ATMOSPHERE MODEL ERROR: Derivatives not defined for atmoshere type %s \n',obj.model)
            end
        end
        
        function plot(obj,fig,lambda,delta)
            %PLOT Verifies atmospehre model
            h = 0:10:11000;
            
            if nargin<=2
                lambda = nan;
                delta  = nan;
            end
            
            if ~exist('fig','var') || isempty(fig)
                figure
            else
                figure(fig)
            end
            
            [rho , P, T , asound] = get_properties(obj,h,lambda,delta);
            
            subplot 221
            
            plot(rho,h/1000,'linewidth',2); hold on
            xlabel('Density [kg/m\textsuperscript{3}]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex')
            grid on
            
            subplot 222
            plot(asound,h/1000,'linewidth',2); hold on
            xlabel('Speed of Sound [m/s]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
            
            subplot 223
            plot(P/1000,h/1000,'linewidth',2); hold on
            xlabel('Pressure [kPa]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
            
            subplot 224
            plot(T,h/1000,'linewidth',2); hold on
            xlabel('Tmperature [K]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
        end
        
        function plotDerivatives(obj,fig)
            %PLOT Verifies atmospehre model
            h = 0:10:11000;
            
            
            if ~exist('fig','var') || isempty(fig)
                figure
            else
                figure(fig)
            end
            
            [drhodh , dpdh, dTdh , dadh] = get_derivatives(obj,h);
            
            subplot 221
            
            plot(drhodh,h/1000,'linewidth',2); hold on
            xlabel('Density Derivative [kg/m\textsuperscript{4}]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex')
            grid on
            
            subplot 222
            plot(dadh,h/1000,'linewidth',2); hold on
            xlabel('Speed of Sound Derivaive [1/s]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
            
            subplot 223
            plot(dpdh/1000,h/1000,'linewidth',2); hold on
            xlabel('Pressure Derivaive [kPa/m]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
            
            subplot 224
            plot(dTdh*1000,h/1000,'linewidth',2); hold on
            xlabel('Tmperature Derivaive [K/km]','interpreter','latex')
            ylabel('Altitude [km]','interpreter','latex') 
            grid on
        end
        
      function out = get.constantTemperature(obj)
		% Computes dependent temperature variable assuming ideal gas law
         out = obj.Hs*obj.g0/obj.gasConstant;
      end
      function out = get.constantPressure(obj)
		% Computes dependent pressure variable assuming ideal gas law
         out = obj.rho0.*obj.constantTemperature*obj.gasConstant;
      end
    end
end

