function [outputy, outputx] = collapseCoordinates(m)
%COLLAPSECOORDINATES returns 2 ordered arrays of the coordinate for each
%white pixel in a 2d binary matrix m
sz=size(m);

outputx=[];
outputy=[];

for y=1:sz(1)
    for x=2:sz(2)
        if m(y,x)==1
            outputx=[outputx x];
            outputy=[outputy y];
        end
    end
end
end

