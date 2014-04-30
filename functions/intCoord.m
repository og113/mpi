%a function which gives the integer lattice coordinates of a point
%argumens are locNum, direction and Nt, ouput is scalar
function XintCoords = intCoords(locNum, direction, Nt)
    global d N;
    intCoordsVec = zeros(d,1);
    Param = zeros(d,1);
    Param(d) = locNum;
    intCoordsVec(d) = floor(locNum/N^(d-2)/Nt);
    for l=1:(d-1-direction)
        Param(d-l) = Param(d-l+1) - intCoordsVec(d-l+1)*N^(d-l-1)*Nt;
        if l~=(d-1)
            intCoordsVec(d-l) = floor(Param(d-l)/N^(d-l-1)/Nt);
        else
            intCoordsVec(d-l) = floor(Param(d-l));
        end
    end
    XintCoords = intCoordsVec(direction+1);
end