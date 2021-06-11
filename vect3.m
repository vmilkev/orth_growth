classdef vect3
    
    properties (GetAccess = private)
        vect_end
        vect_start
    end
    methods
        function plot( obj, r0, r1, col )
            obj.vect_start = r0;
            obj.vect_end = r1;
            q = quiver3(obj.vect_start(1,1),obj.vect_start(1,2),obj.vect_start(1,3),obj.vect_end(1,1),obj.vect_end(1,2),obj.vect_end(1,3),'LineWidth', 3);
            q.Color = col;
        end
    end
    
end