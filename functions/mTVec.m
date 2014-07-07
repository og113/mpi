%gives time in minkowski space as a vector
%argument is Ntm
function XmTVec = mTVec(XNtm)
    global Lb b N;
    xMdim = N*XNtm;
    XmTVec = complex(zeros(xMdim,1));
    for j=0:(xMdim-1)
        XmTVec(j+1) = 1i*Lb - b*intCoord(j,0,XNtm);
    end
end