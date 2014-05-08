%a function which gives the integer lattice coordinates of a point
%only works in 2d
%argumens are locNum, direction and Nt, ouput is scalar
function XintCoord = intCoord(locNum, direction, xNt)
    global N;
    x = floor(locNum/xNt);
    if direction==1
        XintCoord = x;
    else
        XintCoord = locNum - x*xNt;
end