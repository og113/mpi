%a function which gives the integer lattice coordinates of a point
%argumens are locNum and Nt, output is d-vector
function XintCoords = intCoords(locNum, Nt)
    global d N;
    XintCoords = zeros(d,1);
    Param = zeros(d,1);
    Param(d) = locNum;
    XintCoords(d) = floor(locNum/N^(d-2)/Nt);
    for l=1:(d-1)
        Param(d-l) = Param(d-l+1) - XintCoords(d-l+1)*N^(d-l-1)*Nt;
        if l~=(d-1)
            XintCoords(d-l) = floor(Param(d-l)/N^(d-l-1)/Nt);
        else
            XintCoords(d-l) = floor(Param(d-l));
        end
    end
end