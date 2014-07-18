%real t for bubble
%arguments are Edim and Nt
function XrtVec = rtVec(XEdim,XNt)
    global t Lb b;
    XrtVec = zeros(XEdim,1);
    for j=0:(XEdim-1)
        t = intCoord(j,0,XNt);
        XrtVec(j+1) = Lb - b*t;
    end
end