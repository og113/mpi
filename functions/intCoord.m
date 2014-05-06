%a function which gives the integer lattice coordinates of a point
%only works in 2d
%argumens are locNum, direction and Nt, ouput is scalar
function XintCoord = intCoord(locNum, direction, xNt)
    x = floor(locNum/xNt);
    if direction==0
        XintCoord = locNum - xNt*x;
    else
        XintCoord = x;
    end
end