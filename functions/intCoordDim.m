%a function which gives the integer lattice coordinates of a point
%works in arbitrary dimensions
%argumens are locNum, direction and Nt, ouput is scalar
function XintCoordDim = intCoordDim(locNum, direction, xNt)
    global d N;
    intCoordVec = zeros(d,1);
    Param = zeros(d,1);
    Param(d) = locNum;
    intCoordVec(d) = floor(locNum/N^(d-2)/xNt);
    for l=1:(d-1-direction)
        Param(d-l) = Param(d-l+1) - intCoordVec(d-l+1)*N^(d-l-1)*xNt;
        if l~=(d-1)
            intCoordVec(d-l) = floor(Param(d-l)/N^(d-l-1)/xNt);
        else
            intCoordVec(d-l) = floor(Param(d-l));
        end
    end
    XintCoordDim = intCoordVec(direction+1);
end