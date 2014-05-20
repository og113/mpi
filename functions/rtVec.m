%real t for bubble
function XrtVec = rtVec
    global Nt Lt b Edim;
    XrtVec = zeros(Edim,1);
    for j=0:(Edim-1)
        t = intCoord(j,0,Nt);
        XrtVec(j+1) = Lt - b*t;
    end
end