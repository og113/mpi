%a function which gives the integer lattice coordinates of a point
%argumens are locNum, direction and Nt, ouput is scalar
function XintCoords = intCoords(locNum, direction, Nt)
    global d N;
    intCoordsVec = zeros(d,1);
    Param = zeros(d,1);
    Param(d) = locNum;
    intCoordsVec(d) = floor(locNum/N^(d-2)/Nt);
    for k=1:(d-1-direction)
        Param(d-k) = Param(d-k+1) - intCoordsVec(d-k+1)*N^(d-k-1)*Nt;
        if k~=(d-1)
            intCoordsVec(d-k) = floor(Param(d-k)/N^(d-k-1)/Nt);
        else
            intCoordsVec(d-k) = floor(Param(d-k));
        end
    end
    XintCoords = intCoordsVec(direction+1);
end