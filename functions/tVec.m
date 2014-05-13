%like xVec(Nt+Ntm) - gives whole t line - both euclidean and minkowski
function XtVec = tVec;
    global Nt Ntm Lt b Tdim;
    XtVec = zeros(Tdim,1);
    for j=0:(Tdim-1)
        t = intCoord(j,0,Nt+Ntm);
        x = intCoord(j,1,Nt+Ntm);
        if t<Ntm
            XtVec(j+1) = 1i*Lt - b*t;
        else
            XtVec(j+1) = 1i*(Lt - b*(t-Ntm));
        end
    end
end