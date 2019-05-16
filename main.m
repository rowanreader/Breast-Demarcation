function main(file)
pts = inputPts(file);
result = vectors(pts);
out = fopen('output.txt','w');
[n,~] = size(result);
for i=1:1:n
    fprintf(out, '%f %f %f\n', result(i,1), result(i,2), result(i,3));
end
fclose(out);