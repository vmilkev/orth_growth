classdef plane
    
    properties
        center % vectoe
        scale  % scalar
        point1 % vector
        point2 % vector
    end
    
    properties (GetAccess = private)
        precision;
    end
    
    properties (Dependent)
        % structure: coordinates of plane surface
        S
    end
    
    methods

        function S = get.S( obj )
            
            % alocate space for plane
            S.yp = zeros(1,1);
            S.zp = zeros(1,1);
            S.xp = zeros(1,1);
            
            % plane equation:
            % a * x + b * y + c * z = rhs
            
            r1 = obj.point1;
            r2 = obj.point2;

            % vectors on a plane:
            v1 = r1 - obj.center;
            v2 = r2 - obj.center;
                        
            t_range = obj.scale/20;
            s_range = obj.scale/10;
            obj.precision = t_range/100;
            
            i = 1;
            for t = -t_range:obj.precision:t_range
                for s = -s_range:obj.precision:s_range
                    S.xp(i,1) = obj.center(1,1) + v1(1,1) * t + v2(1,1) * s;
                    S.yp(i,1) = obj.center(1,2) + v1(1,2) * t + v2(1,2) * s;
                    S.zp(i,1) = obj.center(1,3) + v1(1,3) * t + v2(1,3) * s;
                    i = i + 1;
                end
            end
        end
        
        function plot( obj, r0, p1, p2, scale )
            obj.center = r0;
            obj.scale = scale;
            obj.point1 = p1;
            obj.point2 = p2;
            
            s = obj.S;
            scatter3( s.xp,s.yp,s.zp, '.' );
        end
    end
    
end