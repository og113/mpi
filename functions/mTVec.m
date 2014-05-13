%gives time in minkowski space as a vector
%argument is Ntm
function XmTVec = mTVec(XNtm)
    global Lt b Mdim;
    XmTVec = zeros(Mdim,1);
    for j=0:(Mdim-1)
        XmTVec(j+1) = 1i*Lt - b*intCoord(j,0,XNtm);
    end
end