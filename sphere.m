classdef sphere
    
    properties
        center
        radius
    end
    
    properties (GetAccess = private)
        precision = 200;
    end
    
    properties (Dependent)
        % structure: coordinates of sphere surface
        S
    end
    
    methods
        %function center = set
        function S = get.S( obj )
            S.x0 = zeros(1,1);
            S.y0 = zeros(1,1);
            S.z0 = zeros(1,1);
            
            phi = 0.0;
            i = 1;
            while phi <= 2*pi
                theta = -pi/2;
                while theta <= pi/2
                    S.x0(i,1) = obj.center(1,1) + obj.radius * cos(theta) * cos(phi);
                    S.y0(i,1) = obj.center(1,2) + obj.radius * cos(theta) * sin(phi);
                    S.z0(i,1) = obj.center(1,3) + obj.radius * sin(theta);
                    theta = theta + pi/obj.precision;
                    i = i + 1;
                end
                phi = phi + pi/obj.precision;
            end
        end
        
        function plot( obj, r0, R )
            obj.center = r0;
            obj.radius = R;
            s = obj.S;
            scatter3( s.x0,s.y0,s.z0, '.' );
        end
    end
    
end