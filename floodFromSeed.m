function floodmask = floodFromSeed(m, point, thresh, floodmask)
%FLOODFROMSEED starting from a seed point, it expands the region until it
%reaches a threshold distance to the border
% input.
% - m: the matrix (binary where 1 represents the border)
% - point: [y,x] array with the coordinates
% - thresh: the threshold distance to the border
% - floodmask: do not assign!!

% floodmask used in recursion, both to keep track of the already explored
% places and to check if we are recursing
if nargin<4
    floodmask=zeros(size(m));
    
    % distance from border 
    m=bwdist(m);
end

%starting from point, expand
if m(point(1), point(2))>thresh
    floodmask(point(1),point(2))=1;
else
    floodmask(point(1),point(2))=-1;
    return
end

% imshow(floodmask) % for debug reasons

trans=[-1,0;0,1;1,0;0,-1];

for i=1:4
    coordinate=point+trans(i,:);

    if coordinate(1) < 1 || coordinate(1) > size(m,1) || coordinate(2) < 1 || coordinate(2) > size(m,2)
        continue
    end

    if floodmask(coordinate(1),coordinate(2))==0
        floodmask=floodmask | floodFromSeed(m, coordinate, thresh, floodmask);
    end

end

