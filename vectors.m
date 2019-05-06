function vectors
t = 10; % degrees
pts = generatePoints;
[center, normal, radius] = CircFit3D(pts);

% find angles between normal and x, and normal and y (degrees)
%theta = atand(normal(3)/normal(1)); % angle between normal and x axis
%phi = atand(normal(3)/normal(2)); % angle between normal and y axis
xVec = [1,0,0];
yVec = [0,1,0];
theta = acos(dot(normal,xVec))/(norm(normal)*norm(xVec));
theta = pi/2 - theta; % adjust for normal
phi = acos(dot(normal,yVec))/(norm(normal)*norm(yVec));
phi = pi/2 - phi;
disp("theta = " + theta);
disp("theta d = " +theta*pi/180);
disp("phi = " + phi);
disp("phi d = "+ phi*pi/180);
% rotate

figure;
alpha = 0:0.01:2*pi;
v=null(normal');
points=repmat(center',1,size(alpha,2))+radius*(v(:,1)*cos(alpha)+v(:,2)*sin(alpha));
plot3(points(1,:),points(2,:),points(3,:),'r-');
hold on;
plot3(pts(:,1), pts(:,2), pts(:,3), 'm.');
plot3([0,10],[0,0],[0,0],'b-');

pts = pts - center;

plot3(pts(:,1), pts(:,2), pts(:,3), 'k.');
plot3([-10*normal(1);10*normal(1)],[-10*normal(2);10*normal(2)],[-10*normal(3);10*normal(3)],'m-');
rotateX = [1 0 0;
          1 cos(theta) sin(theta)
          0 -sin(theta) cos(theta)];
rotateY = [cos(phi) 0 -sin(phi);
          0 1 0
          sin(phi) 0 cos(phi)];
alpha = phi;
rotateZ = [ cos(alpha) sin(alpha) 0;
          -sin(alpha) cos(alpha) 0;
          0 0 1];
pts = pts*rotateX;
%plot3(pts(:,1), pts(:,2), pts(:,3), 'k.');
pts = pts*rotateY;
%figure;
plot3(pts(:,1), pts(:,2), pts(:,3), 'g.');
xlabel("X");
ylabel("Y");
axis square;
%hold on;

% find normal vectors every t degrees
% assume circle is correct and currently on xy plane%%%%%%%%%%
num = 360/t;
% we have center, normal, radius
dist = 20; % how many dots/line
lines = zeros(num*dist,3);
temp = zeros(dist,3);
for i = 0:num
    degree = i*t; % every t degrees
    % make line
    dir = radius/cos(degree);
    unit = dir/norm(dir);
    for j = 1:dist
        temp(j,:) = unit*j;
    end
    lines(i+1:i+dist,:) = temp;
end


% unrotate

% plot

