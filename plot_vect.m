% first = [-0.0169 0.2023 0.9791];
% second = [0.6371 -0.7526 0.1667];
% third = [0.7706 0.1667 -0.1163];
% origin = [0 0 0];
[x1,x2,x3] = sphere;
scale = 0.1;

cameraPos = [0 0 3];
cameraTarget = [0 1 1];
cameraFront = [0 0 -1];
up = [0 1 0];

cameraDirection = (cameraPos - cameraTarget);
cameraDirection = cameraDirection/norm(cameraDirection);
cameraRight = cross(up, cameraDirection);

view = LookAt(cameraPos, cameraTarget, up);
view2 = LookAt(cameraPos, cameraPos-cameraFront, up);

tmp = cameraPos-cameraFront;



figure(1);
box on;
grid on;
hold on;
surf(x1.*scale, x2.*scale, x3.*scale, 'FaceColor', 'y'); % cameraTarget
%surf( (x1+cameraPos(:,1)).*scale, (x2+cameraPos(:,2)).*scale, (x3+cameraPos(:,3)).*scale, 'FaceColor', 'y'); % cameraPos
quiver3( cameraPos(:,1), cameraPos(:,2), cameraPos(:,3), cameraDirection(:,1), cameraDirection(:,2), cameraDirection(:,3), 'Color', 'k' );
quiver3( cameraPos(:,1), cameraPos(:,2), cameraPos(:,3), tmp(:,1), tmp(:,2), tmp(:,3), 'Color', 'r' );
hold off;


% figure(2);
% box on;
% grid on;
% hold on;
% quiver3( origin(:,1), origin(:,2), origin(:,3), first(:,1), first(:,2), first(:,3), 'Color', 'k' );
% quiver3( origin(:,1), origin(:,2), origin(:,3), second(:,1), second(:,2), second(:,3), 'Color', 'r' );
% quiver3( origin(:,1), origin(:,2), origin(:,3), third(:,1), third(:,2), third(:,3), 'Color', 'b' );
% surf(x1.*scale, x2.*scale, x3.*scale, 'FaceColor', 'y');
% hold off;


function viewmatr = LookAt(cameraPos, cameraTarget, up)


cameraDirection = (cameraPos - cameraTarget);
cameraDirection = cameraDirection/norm(cameraDirection);
cameraRight = cross(up, cameraDirection);

matr1 = zeros(4);
matr1(4,4) = 1;
matr1(1, 1) = cameraRight(1,1);
matr1(1, 2) = cameraRight(1,2);
matr1(1, 3) = cameraRight(1,3);
matr1(2, 1) = up(1,1);
matr1(2, 2) = up(1,2);
matr1(2, 3) = up(1,3);
matr1(3, 1) = cameraDirection(1,1);
matr1(3, 2) = cameraDirection(1,2);
matr1(3, 3) = cameraDirection(1,3);
matr2 = eye(4);
matr2(1,4) = -cameraPos(1,1);
matr2(2,4) = -cameraPos(1,2);
matr2(3,4) = -cameraPos(1,3);

viewmatr = matr1*matr2;

end
