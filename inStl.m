function inStl(file)
[v,f,n,name] = stlRead(file);
figure;
plot3(v(:,1),v(:,2),v(:,3),'.b');
axis equal;





%[face,vert,n] = stlread("S1.stl");
%figure;

%[face,vert,n] = stlread("S1.stl");
%patch(face, 'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none','FaceLighting', 'gouraud', 'AmbientStrength', 0.15);
% a = figure;
% face = stlread("S1.stl");
% patch(face(:,1),face(:,2), face(:,3),a);
% axis equal;
% view(30,60);
%,'FaceColor', [0.8,0.8,1], 'EdgeColor', 'none','FaceLighting', 'gouraud', 'AmbientStrength', 0.15);
%camlight('headlight');
%material('dull');
%axis('image');
%view([-135 35]);

% model = createpde;
% importGeometry(model,"S1.stl");
% pdegplot(model,'Facelabels','on')