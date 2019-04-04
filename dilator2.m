function [ dilatedvoxels ] = dilator2(voxels, sz)
% Faster version of dilator
    
    [x, y] = ind2sub(sz,voxels');
    
    xs = [ x-1, x, x+1, x-1,  x+1, x-1, x, x + 1];
    ys = [ y-1, y-1, y-1, y,  y, y+1, y+1, y+1];
    
    inx = (xs > 0) & (xs <= sz(1));
    iny = (ys > 0) & (ys <= sz(2));
    
    in = inx & iny;
    
    dilatedvoxels = unique([voxels', sub2ind(sz, xs(in), ys(in))])';

end