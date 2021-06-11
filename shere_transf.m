
% position of center of sphere
r0 = [1 2 -3];

% radius of sphere
R = 1;

% -------------------------------------------------------------------
% ----------------- construct sphere --------------------------------

% alocate space for sphere surface
x = zeros(1,1);
y = zeros(1,1);
z = zeros(1,1);

precision = 100;
phi = 0.0;
i = 1;
while phi <= 2*pi
    theta = -pi/2;
    while theta <= pi/2
        x(i,1) = r0(1,1) + R * cos(theta) * cos(phi);
        y(i,1) = r0(1,2) + R * cos(theta) * sin(phi);
        z(i,1) = r0(1,3) + R * sin(theta);
        theta = theta + pi/precision;
        i = i + 1;
    end
    phi = phi + pi/precision;
end
% -------------------------------------------------------------------
% -------------------------------------------------------------------


% ----------------- construct plane ---------------------------------
% -------------------------------------------------------------------

% alocate space for plane
yp = zeros(1,1);
zp = zeros(1,1);
xp = zeros(1,1);

i = 1;
for t = -1:0.1:1
    for s = -1:0.1:1
        xp(i,1) = 1-t-s;
        yp(i,1) = t+2*s;
        zp(i,1) = t+2*s;
        i = i + 1;
    end 
end
% -------------------------------------------------------------------
% -------------------------------------------------------------------

figure(1);
box on;
axis equal;
hold on;
scatter3(x,y,z, '.');
scatter3(xp,yp,zp, '.');
xlabel('X');
ylabel('Y');
zlabel('Z');
