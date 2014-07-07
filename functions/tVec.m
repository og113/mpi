%like xVec(Nt+Ntm) - gives whole t line - both euclidean and minkowski
%arguments are Nt, NT and N
function XtVec = tVec(xN,xNa,xNb,xNc)
    global La Lb b;
    xNT = xNa + xNb + xNc;
    XtVec = zeros(xNT*xN,1);
    for j=0:(xNT*xN-1)
        t = intCoord(j,0,xNT);
        x = intCoord(j,1,xNT);
        if t<xNa
            XtVec(j+1) = 1i*Lb + b*(t-xNa);
        elseif t<(xNa+xNb)
            XtVec(j+1) = 1i*(Lb - b*(t-xNa));
        else
            XtVec(j+1) = b*(t-xNa-xNb);
        end
    end
end