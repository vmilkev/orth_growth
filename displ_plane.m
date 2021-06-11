clear;
% -------------------------------------------------------------------
clf;
figure(1);
box on;
axis equal;
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% position of the center of a first (free) sphere (cluster initiator)
s11_r0 = [0.00001 0.0 0.0];
% radius of the first sphere
s11_R = 3;

eps = 0.5; % shift of new sphere from the common origin

% ------------------------------------------------------------------


% ----------------- FIRST DIVISION ----------------------------------

% random [x y] point that will determine a first division plane
orthdir = [(-1+2*rand)*10 (-1+2*rand)];

% get random normal to a first division plane
%            a          b            c
s1_n = [(-1+2*rand) (-1+2*rand) (-1+2*rand)];

% get vectors to arbitrary points on the plane
% should satisfy equation: a*x + b*y + c*z = rhs.
rhs = dot(s11_r0,s1_n);
s1_r1 = [orthdir (rhs - s1_n(1,1)*orthdir(1,1) - s1_n(1,2)*orthdir(1,2))/s1_n(1,3)]; % vector to arbitrary point 1 in a FIRST DIVISION plane
s1_r2 = [-orthdir (rhs + s1_n(1,1)*orthdir(1,1) + s1_n(1,2)*orthdir(1,2))/s1_n(1,3)]; % vector to arbitrary point 2 in a FIRST DIVISION plane

% normal to the first plane at the point s1_r0 (scaled by R of the initial sphere)
n1 = getNormal( s11_r0, s1_r1, s1_r2 ).*s11_R;
% visualize:
% 1) first sphere
% 2) first division plane
% 3) directions of growth

% s1 = sphere;
% p1 = plane;
% vn1 = vect3;
% s1.plot(s11_r0, s11_R);
%p1.plot(s1_r0, s1_r1, s1_r2, s1_R);
% vn1.plot(s11_r0, n1);
% vn1.plot(s11_r0, -n1);
% -------------------------------------------------------------------

% ----------------- SECOND DIVISION ---------------------------------

% get INITIAL positions for new spheres
s12_r0 =  getNewSphere( s11_r0, n1, eps );%s11_r0 - n1 * eps1;
s13_r0 =  getNewSphere( s11_r0, n1, -eps );%s11_r0 - n1 * eps2;

% position for the second division plane given the centers of the two new
% spheres
p2_pos = getPlanePts2( s12_r0, s13_r0 );
p3_pos = getPlanePts2( s13_r0, s12_r0 );

n12 = getNormal( p2_pos.one, p2_pos.two, p2_pos.three ).*s11_R;
n13 = getNormal( p3_pos.one, p3_pos.two, p3_pos.three ).*s11_R;

% visualize:
% 1) next sphere
% 2) second division plane
% 3) directions of growth

% s12 = sphere;
% s13 = sphere;
% p21 = plane;
% p22 = plane;
% vn12 = vect3;
% vn13 = vect3;
% s12.plot(s12_r0, s11_R);
% s13.plot(s13_r0, s11_R);
% p21.plot(p2_pos.one, p2_pos.two, p2_pos.three, s1_R);
% p22.plot(p3_pos.one, p3_pos.two, p3_pos.three, s1_R);
% vn12.plot(s12_r0, n12);
% vn12.plot(s12_r0, -n12);
% vn13.plot(s13_r0, n13);
% vn13.plot(s13_r0, -n13);

% -------------------------------------------------------------------

% ----------------- THIRD DIVISION ----------------------------------

% division of sphere s2:
% get INITIAL positions for new spheres
s21_r0 = getNewSphere( s12_r0, n12, eps );%s12_r0 - n12 * eps1;
s22_r0 = getNewSphere( s12_r0, n12, -eps );%s12_r0 - n12 * eps2;

% s21 = sphere;
% s22 = sphere;
% s21.plot(s21_r0, s11_R);
% s22.plot(s22_r0, s11_R);

% division of sphere s3:
% get INITIAL positions for new spheres
s31_r0 = getNewSphere( s13_r0, n13, eps );%s13_r0 - n13 * eps1;
s32_r0 = getNewSphere( s13_r0, n13, -eps );%s13_r0 - n13 * eps2;

% s31 = sphere;
% s32 = sphere;
% s31.plot(s31_r0, s11_R);
% s32.plot(s32_r0, s11_R);

% new division directions

% 1) division of S21 (the point of growth):
growthfrom = s21_r0;
contact1 = s31_r0;
contact2 = s22_r0;
scaler = s11_R;
n21 = getCross( contact1-growthfrom, contact2-growthfrom ).*scaler;
vn21 = vect3;
% vn21.plot( growthfrom, n21 );
% vn21.plot( growthfrom, -n21 );
% sphere S21 disappears and two new S211 & S212 appears
% get INITIAL positions for new spheres
s211_r0 = getNewSphere( growthfrom, n21, eps );%s21_r0 - n21 * eps1;
s212_r0 = getNewSphere( growthfrom, n21, -eps );%s21_r0 - n21 * eps2;


% 2) division of S22:
growthfrom = s22_r0;
contact1 = s21_r0;
contact2 = s32_r0;
scaler = s11_R;
n22 = getCross( contact1-growthfrom, contact2-growthfrom ).*scaler;
vn22 = vect3;
% vn22.plot( growthfrom, n22 );
% vn22.plot( growthfrom, -n22 );
% sphere S22 disappears and two new S221 & S222 appears
% get INITIAL positions for new spheres
s221_r0 = getNewSphere( growthfrom, n22, eps );%s22_r0 - n22 * eps1;
s222_r0 = getNewSphere( growthfrom, n22, -eps );%s22_r0 - n22 * eps2;


% 3) division of S31:
growthfrom = s31_r0;
contact1 = s21_r0;
contact2 = s32_r0;
scaler = s11_R;
n31 = getCross( contact1-growthfrom, contact2-growthfrom ).*scaler;
vn31 = vect3;
% vn31.plot( growthfrom, n31 );
% vn31.plot( growthfrom, -n31 );
% sphere S31 disappears and two new S311 & S312 appears
% get INITIAL positions for new spheres
s311_r0 = getNewSphere( growthfrom, n31, eps );%s31_r0 - n31 * eps1;
s312_r0 = getNewSphere( growthfrom, n31, -eps );%s31_r0 - n31 * eps2;


% 4) division of S32:
growthfrom = s32_r0;
contact1 = s22_r0;
contact2 = s31_r0;
scaler = s11_R;
n32 = getCross( contact1-growthfrom, contact2-growthfrom ).*scaler;
vn32 = vect3;
% vn32.plot( growthfrom, n32 );
% vn32.plot( growthfrom, -n32 );
% sphere S32 disappears and two new S321 & S322 appears
% get INITIAL positions for new spheres
s321_r0 = getNewSphere( growthfrom, n32, eps );%s32_r0 - n32 * eps1;
s322_r0 = getNewSphere( growthfrom, n32, -eps );%s32_r0 - n32 * eps2;


sph = sphere;
% sph.plot(s211_r0, 3);
% sph.plot(s212_r0, 3);
% sph.plot(s221_r0, 3);
% sph.plot(s222_r0, 3);
% sph.plot(s311_r0, 3);
% sph.plot(s312_r0, 3);
% sph.plot(s321_r0, 3);
% sph.plot(s322_r0, 3);
%--------------------------------------------------------------------

% ----------------- Test Algorithm ----------------------------------
c = vect3;
% get list (coordinates of centers) of all existing (live) cells
R = 3;
% all_cells(1,:) = s211_r0;
% all_cells(2,:) = s212_r0;
% all_cells(3,:) = s221_r0;
% all_cells(4,:) = s222_r0;
% all_cells(5,:) = s311_r0;
% all_cells(6,:) = s312_r0;
% all_cells(7,:) = s321_r0;
% all_cells(8,:) = s322_r0;

s211_r0 = [0 0 0];
s212_r0 = [3 0 0];
s221_r0 = [6 0 0];
s222_r0 = [3 3 0];
s311_r0 = [3 0 3];
s312_r0 = [3 0 -3];
all_cells(1,:) = s211_r0;
all_cells(2,:) = s212_r0;
all_cells(3,:) = s221_r0;
all_cells(4,:) = s222_r0;
all_cells(5,:) = s311_r0;
all_cells(6,:) = s312_r0;
% sph.plot(s211_r0, 3);
% sph.plot(s212_r0, 3);
% sph.plot(s221_r0, 3);
% sph.plot(s222_r0, 3);
% sph.plot(s311_r0, 3);
% sph.plot(s312_r0, 3);

g1 = mcell(1,s211_r0,R);
g2 = mcell(2,s212_r0,R);
g3 = mcell(3,s221_r0,R);
g4 = mcell(4,s222_r0,R);
g5 = mcell(5,s311_r0,R);
g6 = mcell(6,s312_r0,R);

g7 = mcell(7,[20 20 20],R);

cluster{1,1} = g1;
cluster{2,1} = g2;
cluster{3,1} = g3;
cluster{4,1} = g4;
cluster{5,1} = g5;
cluster{6,1} = g6;
cluster{7,1} = g7;

% g_nov = g7.growth( cluster, 4 );
%g_nov
g_nov = g7.growth( cluster, 5 );
g7.vplot;
% g_nov

%g7.splot;
g7.vplot;
g_nov.splot;

cluster{8,1} = g_nov;

g_n = g_nov.growth(cluster, 5);
%g_n.splot;
g_nov.vplot;
%g_nov.splot;
cellNum = size(all_cells,1);

for i = 1:cellNum
    
    % get contact list for the cell
    cnt = getContacts( i, all_cells, R );
    contactNum = size(cnt, 1);
            
    origin = all_cells(i,:);
    
    if (contactNum == 1)
        % get growth direction
        confinement = cnt(1,:);
        orthdir = getOrthVect( confinement );
        growthdir = getCross(confinement,orthdir);
        
        newcell_1 = getNewSphere( origin,growthdir, 0.0 );
        newcell_2 = getNewSphere( origin,growthdir, 0.5 );
        
        % check growth direction vector
        c.plot(origin,growthdir.*2.5*R);
%         c.plot(origin,confinement.*1.5*R);
%         dot(growthdir,confinement)
    end
    
%     if (contactNum == 20)
%         
%         confinement1 = cnt(1,:);
%         confinement2 = cnt(2,:);
%         
%         % check if the confinement vectors are not parallel
%         if ( abs( dot(confinement1,confinement2) ) >= 1.0 )
%             orthdir = getOrthVect( confinement1 );
%             growthdir = getCross(confinement1,orthdir);
%             newcell_1 = getNewSphere( origin,growthdir, 0.5 );
%             newcell_2 = getNewSphere( origin,growthdir, -0.5 );
% 
%             % check growth direction vector
%             c.plot(origin,growthdir.*2.5*R);
% %             c.plot(origin,confinement1.*1.5*R);
% %             c.plot(origin,confinement2.*1.5*R);
%             dot(growthdir,confinement1)
%         else
%             growthdir = getCross(confinement1,confinement2);
%             newcell_1 = getNewSphere( origin,growthdir, 0.5 );
%             newcell_2 = getNewSphere( origin,growthdir, -0.5 );
% 
%             % check growth direction vector
% %             c.plot(origin,confinement1.*1.5*R);
% %             c.plot(origin,confinement2.*1.5*R);
%             c.plot(origin,growthdir.*2.5*R);
%             
%             dot(growthdir,confinement1)
%             dot(growthdir,confinement2)
% 
%         end
%     end
    
    if (contactNum >= 2)
        for i = 1:contactNum
            %c.plot(origin,cnt(i,:).*1.0*R);
        end
        % loop over contacts and find all possible growth directions
        l = 1;
        for i = 1:contactNum-1
            for j = i+1:contactNum
                confinement1 = cnt(i,:);
                confinement2 = cnt(j,:);
                dir = getCross(confinement1,confinement2);
                if ( norm(dir) ~= 0.0 ) % if norm(dir) == 0.0 than conf1 is parallel to conf2
                    tmp_directions(l,:) = dir;
                    tmp_directions(l+1,:) = -dir;
%                     c.plot(origin,dir.*1.0*R);
%                     c.plot(origin,-dir.*1.0*R);
                    l = l + 2;
                end
            end
        end
        % loop over directions and remove those directed towards
        % existing confinements
        if ~exist('tmp_directions','var')            
            orthdir = getOrthVect( cnt(1,:) );
            tmp_directions(1,:) = getCross(cnt(1,:),orthdir);
        end
        
        dirNum = size(tmp_directions,1);
        l = 1;
        for i = 1:dirNum
            dir = tmp_directions(i,:);
            accepted = true;
            for j = 1:contactNum
                conf = cnt(j,:);
                check = dot(dir,conf);
                if ( check < -0.9 ) % if check == -1 growth directed towards confinement
                    %c.plot(origin,dir.*1.0*R);
                    accepted = false;
                    continue;
                end
            end
            if ( accepted )
                directions(l,:) = dir;
                l = l + 1;
            end
        end
        
        % allow only one direction randomly
        dirNum = size(directions,1);
        growthdir = directions( randi(dirNum),: );
        c.plot(origin,growthdir.*2.5*R);
    end

end

%------------------------------------------------------------


function n = getNormal( p0, p1, p2 )

% calculate normale to a plane at the point p0
% using two other (arbitrary) points on a plane: p1 and p2

% vectors on a plane:
v1 = p1 - p0; % vector laying on a plane from p0 to p1
v2 = p2 - p0; % vector laying on a plane from p0 to p2

normale = [ (v1(1,2)*v2(1,3) - v1(1,3)*v2(1,2)) -(v1(1,1)*v2(1,3) - v1(1,3)*v2(1,1)) (v1(1,1)*v2(1,2) - v1(1,2)*v2(1,1)) ];

% make a unit vector
n = normale./norm(normale);

end

function v = getCross( v1, v2 )

v1 = [ (v1(1,2)*v2(1,3) - v1(1,3)*v2(1,2)) -(v1(1,1)*v2(1,3) - v1(1,3)*v2(1,1)) (v1(1,1)*v2(1,2) - v1(1,2)*v2(1,1)) ];

% make a UNIT vector
if ( norm(v1) == 0.0 )
    v = v1;
else
    v = v1./norm(v1);
end

end

function r2 = getOrthVect( r1 )

if ( r1(1,1) == 0.0 )
    r1(1,1) = 0.00001;
end

if ( r1(1,2) == 0.0 )
    r1(1,2) = 0.00001;
end

r2 = [1/r1(1,1) -1/r1(1,2) 0.0];

end

function p = getPlanePts2( r1, r2 )

%     returns three points on a plane
%     the first two points are known (probably centers of two spheres)
%     and the third one is an arbitrary point (determines plane rotation along r1->r2 vector)

p.one = r1;
p.two = r2;
x = (-1+2*rand);
y = (-1+2*rand);
z = (-1+2*rand);
p.three = [ x y z ];
end


function r = getNewSphere( origin, direction, shift)
r = origin - direction * shift;
end

function contact = getContacts( whichcell, cells, radius )

% get contact list

cellNum = size(cells,1);

i = whichcell;
    l = 1;
    for j = 1:cellNum
        if ( i ~= j )
            v1 = cells(i,:);
            v2 = cells(j,:);
            v12 = v1-v2;
            d = norm(v12);
            if ( (radius + radius - d) >= 0.7*(radius + radius)/2 )
                contact(l,:) = v12./norm(v12);
                l = l + 1;
            end
        end
    end
end
