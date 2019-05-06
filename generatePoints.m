% function to generate n points in the shape of a circle or radius r randomly in space
% just a test function to generate input
function pts = generatePoints
n = 100; % number of points
r = 50; % radius
x0 = 0; % Center in x
y0 = 0; % Center in y
z0 = 0; % Center in z 
error = 5; % max error
transCoord = 100; % range for transformation starting point
% create the set of points
pts = zeros(n,3); % n rows, 3 columns

alpha = zeros(n,1);
for i = 1:n
    % get random angle
    angle = rand(1)*2*pi;

    % find [x,y,z] of point in that angle of radius r
    pt = [cos(angle)*r + x0, sin(angle)*r + y0, 0 + z0];
    % add error to all [x,y,z]
    xErr = randi([-error,error],1,1);
    yErr = randi([-error,error],1,1);
    zErr = randi([-error,error],1,1);
    pt = [pt(1) + xErr, pt(2) + yErr, pt(3) + zErr];
    pts(i,:) = pt;    
    alpha(i,:) = angle;
end

% original circle
%plot3(pts(:,1), pts(:,2), pts(:,3), 'bo');
%hold on;

% transform circle to random rotation about x axis and y axis
% random angles
theta = randi([0,360], 1, 1);
theta_rad = pi() * theta / 180.0;
phi = randi([0,360], 1, 1);
phi_rad = pi() * phi / 180.0;

% rotation matrices
rotateX = [
    1 0 0;
    0 cos(theta_rad) -sin(theta_rad);
    0 sin(theta_rad) cos(theta_rad)];
rotateY = [
    cos(phi_rad) 0 sin(phi_rad);
    0 1 0;
    -sin(phi_rad) 0 cos(phi_rad)];

for i = 1: n
    pts(i,:) = pts(i,:) *rotateX';
    pts(i,:) = pts(i,:) *rotateY';
end

% transform circle to random 3D space
spot = randi([-transCoord, transCoord],1,3);
for i = 1: n
    pts(i,:) = pts(i,:) + spot;
end

% figure;
% plot3(pts(:,1), pts(:,2), pts(:,3), 'r.');
% axis equal;
% hold on;
% 
% [center, normal, rad] = CircFit3D(pts);
% alpha = 0:0.01:2*pi;
% v=null(normal');
% points=repmat(center',1,size(alpha,2))+rad*(v(:,1)*cos(alpha)+v(:,2)*sin(alpha));
% plot3(points(1,:),points(2,:),points(3,:),'r-');
%axis square;
end
