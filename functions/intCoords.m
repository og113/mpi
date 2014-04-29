%a function which gives the integer lattice coordinates of a point
%argumens are locNum and Nt, output is d-vector
function XintCoords = intCoords(locNum, Nt)
    global d N;
    XintCoords = zeros(d,1);
    Param = zeros(d,1);
    Param(d) = locNum;
    XintCoords(d) = floor(locNum/N^(d-2)/Nt);
    for k=1:(d-1)
        Param(d-k) = Param(d-k+1) - XintCoords(d-k+1)*N^(d-k-1)*Nt;
        if k~=(d-1)
            XintCoords(d-k) = floor(Param(d-k)/N^(d-k-1)/Nt);
        else
            XintCoords(d-k) = floor(Param(d-k));
        end
    end
end