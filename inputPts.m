function pts = inputPts(file)

fileID = fopen(file,'r');
a = fscanf(fileID, '%f');
fclose(fileID);
[n,~] = size(a);
m = n/3; % actual number of points
pts = zeros(m, 3);
for i=0:m-1
    pts(i+1,:) = [a(3*i+1),a(3*i+2),a(3*i+3)];
end

end