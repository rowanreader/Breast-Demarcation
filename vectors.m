function unrotated = vectors(pts)
t = 10; % degrees of separation
% pts = generatePoints;
[center, normal, radius] = CircFit3D(pts);
plot = 1;
% find angles between normal and x, and normal and y (degrees)
%theta = atand(normal(3)/normal(1)); % angle between normal and x axis
%phi = atand(normal(3)/normal(2)); % angle between normal and y axis

zVec = [0,0,1];
% get angle between normal and zVec
angle = acos(dot(zVec, normal)/(norm(zVec)*norm(normal)));

% draw circle
if plot == 1
    figure;
    alpha = 0:0.01:2*pi;
    v=null(normal');
    points=repmat(center',1,size(alpha,2))+radius*(v(:,1)*cos(alpha)+v(:,2)*sin(alpha));
    plot3(points(1,:),points(2,:),points(3,:),'k-');
    hold on;
    plot3(pts(:,1), pts(:,2), pts(:,3), 'm.');
end
% translate points to center
pts = pts - center;
normal = normal/(norm(normal));
% rotate about x, then y => also make 'unrotate' matrices
a = normal(1);
b = normal(2);
c = normal(3);
dir = sqrt(normal(2).^2 + normal(3).^2);
l = sqrt(a^2 + b^2 + c^2);
v = sqrt(b^2 + c^2);
cos_t = c/v;
sin_t = b/v;
rotate_x = [1,0,0;...
    0,cos_t,-sin_t;...
    0, sin_t, cos_t];
rotate_xi = [1,0,0;...
    0,cos_t, sin_t;...
    0, -sin_t, cos_t];
rotate_y = [dir, 0, -a;...
    0,1,0;...
    a, 0, dir];
rotate_yi = [dir, 0, a;...
    0,1,0;...
    -a, 0, dir];

% in theory, rotate around z too, but it's a circle... somewhat
% unneccessary
rotate_z = [cosd(angle), -sind(angle), 0;...
    sind(angle), cosd(angle), 0;...
    0,0,1];
pts = (rotate_z*(rotate_y*(rotate_x*pts')));
%figure;
zNorm = (rotate_z*(rotate_y*(rotate_x*normal))); % should be [0,0,1], since we rotated to be parallel to z axis
% draw rotated circle and points
if plot == 1
    v = null(zNorm');
    circle=repmat([0;0;0],1,size(alpha,2))+radius*(v(:,1)*cos(alpha)+v(:,2)*sin(alpha));
    plot3(circle(1,:),circle(2,:),circle(3,:),'r-');
    pts = pts';
    plot3(pts(:,1), pts(:,2), pts(:,3), 'g.');
    xlabel("X");
    ylabel("Y");
end


% find normal vectors every t degrees
num = 360/t;
% we have center, normal, radius
dirNum = 30; % how many dots/line
space = 5; % separation of dots
temp = zeros(dirNum,3); % holds points for a given line
lines = zeros(num*dirNum,3); % holds points for all lines
for i = 0:num
    degree = i*t; % every t degrees
    % find point on circle at degree
    x = radius*sind(degree);
    y = radius*cosd(degree);    
    % find perpindicular to tangent (literally just treat as direction
    % vector)
    vector = [x,y];
    vector = vector/norm(vector); % normalize
    % find dirNum points along line, starting radius away, record
    for j = 1:dirNum
        % since at origin and in line with z axis, all z coordinates are 0
        temp(j,:) = [vector*((j*space)+radius),0];
    end
    lines(i*dirNum+1:(i+1)*dirNum,:) = temp; 
end
% plot on rotated circle
if plot == 1
    plot3(lines(:,1), lines(:,2),lines(:,3),'k.');
    axis equal;
end
% unrotate
unrotated = center' + (rotate_xi*(rotate_yi*lines'));

% plot on original circle
if plot == 1
    plot3(unrotated(1,:), unrotated(2,:),unrotated(3,:),'b.');
end
unrotated = unrotated';
end

