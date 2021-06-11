classdef mcell < handle
    
    properties
        origin
        radius
        name
        growthdir
        isElongating
        growthRate
        % for controlling of growth
        isLinked
        linkedID
        % geometrical properties
        conf_degree
        conf_num
        age
        currentElongation
        confinements
        
        isVirtual
    end
    
    properties (GetAccess = private)
        %age
        %conf_num
        %confinements
        origin_nov
        growthRate_nov
        %currentElongation
        % these are for visualisation
        S
        V
        divisiondir
        maxElongation
        confLimit
    end
    
    properties (Constant, GetAccess = private)
        % age at which cell starts to growth
        initElongationTime = 5.0;
        growthRateConstant = 0.1;
        % relative distance between growing cells
        % should be related to cell's radius
        %maxElongation = 3.50;
        % max value for random ID number
        maxID = 10000;
    end
    
    methods
        function this = mcell( orig, rad )
            rng('shuffle');
            if nargin == 2
                this.origin = orig;
                this.origin_nov = [];
                this.radius = rad;
                this.name = randi(this.maxID);
                this.age = 0.0;
                this.conf_num = 0;
                this.confinements = [];
                this.growthdir = [];
                this.divisiondir = [];
                this.S = sphere;
                this.V = vect3;
                %this.initDivision = false;
                this.isLinked = false;
                this.isElongating = false;
                this.linkedID = [];
                this.currentElongation = 0.0;
                this.growthRate = this.growthRateConstant;
                this.conf_degree = 0.0;
                
                this.isVirtual = false;
                this.maxElongation = 0.7*this.radius;
                this.confLimit = 1.0;
            end
        end
        
        function splot( this )
            this.S.plot( this.origin, this.radius*0.2 );
        end
        function vplot_dir( this, val )
            if ( ~isempty(this.divisiondir) )
                this.V.plot( this.origin, val*this.divisiondir.*this.radius*0.5, 'black' );
            end
        end
        function vplot_conf( this,i_conf )
            if ( ~isempty(this.confinements) )
                % here, '-confinements' is just for convinience, so arrow
                % points towards confinement (in the code it is opposite)
                this.V.plot( this.origin, -this.confinements(i_conf,:).*this.radius*0.5, 'red' );
            end
        end

        function newcell = addcell( this )
            if ( this.isElongating && ~this.isLinked )
                % a new cell
                newcell = mcell(this.origin_nov,this.radius);
                % the new cell is still a part of an original cell
                newcell.isLinked = true;
                newcell.linkedID = this.name;
                newcell.isElongating = true;
                newcell.growthdir = this.growthdir;
                newcell.growthRate = this.growthRate_nov;
                newcell.age = this.age;
                
                newcell.isVirtual = true;
                % change the state of the original cell
                this.isLinked = true;
                this.linkedID = newcell.name;
                this.origin_nov = [];
                
            end
        end
        
        function growth( this, cells, time )
            this.age = this.age + time;
            if ( this.isElongating )
                
                if ( ~this.isLinked )
                    return;
                end
                                
                % check whether to release cells
                %if (this.currentElongation >= this.maxElongation)
                if ( this.currentElongation/this.maxElongation >= 1.0 )
                    % release the cell and reset its properties
                    this.isLinked = false;
                    this.age = 0.0;
                    this.growthdir = [];
                    this.conf_num = 0;
                    this.isElongating = false;
                    %this.linkedID = [];
                    this.currentElongation = 0.0;
                    this.origin_nov = [];
                    this.confinements = [];
                    this.divisiondir = [];
                    this.growthRate = this.growthRateConstant;
                    
                    this.isVirtual = false;
                    return;
                end
                
                % correct growth rate to be proportional
                % to the cell's confinement
                %this.growthRate = this.growthRate;
                % move origin of the cell
                this.origin = this.getNewOrigin( this.origin, this.growthdir, this.growthRate );
                % calculate elongation distance
                cellNum = size(cells,1);
                for i = 1:cellNum
                    if ( this.name == cells{i,1}.linkedID )
                        if ( this.isVirtual )
                            this.currentElongation = cells{i,1}.currentElongation;
                            break;
                        end
                        v1 = this.origin;
                        v2 = cells{i,1}.origin;
                        v12 = v1-v2;
                        this.currentElongation = norm(v12);
                        break;
                    end
                end
            elseif ( this.age >= this.initElongationTime )
                % cell is ready (mature enough) for elongation
                this.getGrowthDirections( cells );
                
                if (this.conf_degree > 0.9)
                    this.conf_degree = 1.0;
                end

                switch this.conf_num
                    case {0, 1, 2}
                        this.divisiondir = this.growthdir;
                        currorigin = this.origin;
                        this.growthRate = this.growthRateConstant * (1.0-this.conf_degree);
                        this.growthRate_nov = -this.growthRateConstant * (1.0-this.conf_degree);
                    otherwise
                        this.divisiondir = this.growthdir;
                        currorigin = this.origin;
                        this.growthRate = 0.0;
                        this.growthRate_nov = -this.growthRateConstant*2.0 * (1.0-this.conf_degree);
                end
                % initial elongation event
                this.origin = this.getNewOrigin( currorigin, this.growthdir, this.growthRate );
                this.origin_nov = this.getNewOrigin( currorigin, this.growthdir, this.growthRate_nov );
                this.isElongating = true;
            end
            
        end
    end
    
    methods ( Access = private )
        
        function n = getNormal( this, p0, p1, p2 )
            if ( norm(p0) == 0.0 )
                p0 = p0 + 0.0001;
            end
            % calculate normale to a plane at the point p0
            % using two other (arbitrary) points on a plane: p1 and p2
            % vectors on a plane:
            v1 = p1 - p0; % vector laying on a plane from p0 to p1
            v2 = p2 - p0; % vector laying on a plane from p0 to p2
            normale = [ (v1(1,2)*v2(1,3) - v1(1,3)*v2(1,2)) -(v1(1,1)*v2(1,3) - v1(1,3)*v2(1,1)) (v1(1,1)*v2(1,2) - v1(1,2)*v2(1,1)) ];
            % make a unit vector
            if ( norm(normale) == 0.0 )
                n = normale;
            else
                n = normale./norm(normale);
            end
        end
        
        function v = getCross( this, v1, v2 )
            v1 = [ (v1(1,2)*v2(1,3) - v1(1,3)*v2(1,2)) -(v1(1,1)*v2(1,3) - v1(1,3)*v2(1,1)) (v1(1,1)*v2(1,2) - v1(1,2)*v2(1,1)) ];
            % make a UNIT vector
            if ( norm(v1) == 0.0 )
                v = v1;
            else
                v = v1./norm(v1);
            end
        end
        
        function r2 = getOrthVect( this, r1 )
            if ( r1(1,1) == 0.0 )
                r1(1,1) = 0.00001;
            end
            if ( r1(1,2) == 0.0 )
                r1(1,2) = 0.00001;
            end
            r2 = [1/r1(1,1) -1/r1(1,2) 0.0];
        end
        
        function getConfinements( this,cells )                            
            cellNum = size(cells,1);
            t_confinements = zeros(1,3);
            modified = false;
            l = 1;
            for j = 1:cellNum
                %if ( (this.name ~= cells{j,1}.name) && ~cells{j,1}.isVirtual )
                if ( (this.name ~= cells{j,1}.name) )
                    v1 = this.origin;
                    v2 = cells{j,1}.origin;
                    v12 = v1-v2;
                    d = norm(v12);
                    %if ( (this.radius + cells{j,1}.radius - d) >= 0.7*(this.radius + cells{j,1}.radius)/2 )
                    if ( d < this.radius*this.confLimit )                        
                        %this.confinements(l,:) = v12./norm(v12);
                        t_confinements(l,:) = v12./norm(v12);
                        modified = true;
                        l = l + 1;
                    end
                end
            end
            
            % select only orthogonal directions
            if (modified)
                if ( size(t_confinements,1) > 2 )
                    l = 1;
                    for ii = 1:size(t_confinements,1)
                        for jj = 1:size(t_confinements,1)
                            if (ii ~= jj)
                                res = abs( dot( t_confinements(ii,:),t_confinements(jj,:) ) );
                                if ( res < 0.05 )
                                    % if the confinement vectors are
                                    % orthogonal:
                                %if ( res > 0.99 || res < 0.05 )
                                    this.confinements(l,:) = t_confinements(ii,:);
                                    l = l + 1;
                                    break;
                                end
                            end
                        end
                    end
                else
                    this.confinements = t_confinements;
                end
            end
            
            this.conf_num = size(this.confinements,1);
            this.conf_degree = 1-1/l;
            %disp([size(this.confinements,1)]);
        end
        
        function getGrowthDirections( this,cells )
            % get contact list for the cell
            this.getConfinements( cells );
            
            if ( this.conf_num == 0 )
                % random [x y] point that will determine a first division plane
                orthdir = [(-1+2*rand)*10 (-1+2*rand)];
                % get random normal to a first division plane
                %            a          b            c
                nvect = [(-1+2*rand) (-1+2*rand) (-1+2*rand)];
                % get vectors to arbitrary points on the plane
                % should satisfy equation: a*x + b*y + c*z = rhs.
                rhs = dot(this.origin,nvect);
                r1 = [orthdir (rhs - nvect(1,1)*orthdir(1,1) - nvect(1,2)*orthdir(1,2))/nvect(1,3)]; % vector to arbitrary point 1 in a FIRST DIVISION plane
                r2 = [-orthdir (rhs + nvect(1,1)*orthdir(1,1) + nvect(1,2)*orthdir(1,2))/nvect(1,3)]; % vector to arbitrary point 2 in a FIRST DIVISION plane
                % normal to the first plane at the point s1_r0 (scaled by R of the initial sphere)
                this.growthdir = this.getNormal( this.origin, r1, r2 );
            elseif (this.conf_num == 1)
                %this.confinements(:,1)
                orthdir = this.getOrthVect( this.confinements );
                this.growthdir = this.getCross(this.confinements,orthdir);
            else
                % loop over contacts and find all possible growth directions
                tmp_directions = [];
                directions = [];
                l = 1;
                for i = 1:this.conf_num-1
                    for j = i+1:this.conf_num
                        confinement1 = this.confinements(i,:);
                        confinement2 = this.confinements(j,:);
                        dir = this.getCross(confinement1,confinement2);
                        if ( norm(dir) ~= 0.0 ) % if norm(dir) == 0.0 than conf1 is parallel to conf2
                            tmp_directions(l,:) = dir;
                            l = l + 1;
%                             tmp_directions(l+1,:) = -dir;
%                             l = l + 2;
                        end
                    end
                end
                % loop over directions and remove those directed towards
                % existing confinements
                if ( isempty(tmp_directions) )
                    orthdir = this.getOrthVect( this.confinements(1,:) );
                    tmp_directions(1,:) = this.getCross(this.confinements(1,:),orthdir);
                end
                
                dirNum = size(tmp_directions,1);
                l = 1;
                for i = 1:dirNum
                    dir = tmp_directions(i,:);
                    accepted = true;
                    for j = 1:this.conf_num
                        check = dot( dir,this.confinements(j,:) );
                        if ( check < -0.9 ) % if check == -1 growth directed towards confinement
                            accepted = false;
                            break;
                        end
                    end
                    if ( accepted )
                        directions(l,:) = dir;
                        l = l + 1;
                    end
                end
                % allow only one direction randomly
                dirNum = size(directions,1);
                this.growthdir = directions( randi(dirNum),: );
            end
        end
        
        function neworigin = getNewOrigin( this, currentorigin, growthdirection, cellshift)
            neworigin = currentorigin - growthdirection * cellshift;
        end
        
    end
    
end