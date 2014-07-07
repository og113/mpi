%like xVec(Nt+Ntm) - gives whole t line - both euclidean and minkowski
%arguments are Nt, NT and N
function XtVec = tVec(xNt,xNT,xN)
    global Lt b Ltm;
    XtVec = zeros(xNT*xN,1);
    for j=0:(xNT*xN-1)
        t = intCoord(j,0,xNT);
        x = intCoord(j,1,xNT);
        if t<(xNT-xNt)
            XtVec(j+1) = 1i*Lt + b*(t-xNT+xNt);
        else
            XtVec(j+1) = 1i*(Lt - b*(t-(xNT-xNt)));
        end
    end
end