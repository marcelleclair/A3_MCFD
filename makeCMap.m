function cMap = makeCMap(l,w,ds,boxes)

nx = ceil(l/ds) + 1;
ny = ceil(w/ds) + 1;
cMap = ones(ny,nx);

x = ds*(0:1:nx-1);
X = repmat(x,[ny 1]);
X = X(:);
y = transpose(ds*(0:1:ny-1));
Y = repmat(y,[1 nx]);
Y = Y(:);
p = transpose([X Y]);

for b = 1 : length(boxes)
    IN = isInside(p,boxes(b));
    IN_mat = reshape(IN,ny,nx);
    cMap(IN_mat) = boxes(b).sig;
end

end