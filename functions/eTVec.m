%gives time in euclidean space as a vector
%arguments are Nt and N
function XeTVec = eTVec(xNb, xN)
    global Lb b;
    XeTVec = zeros(xNb*xN,1);
    for j=0:(xNb*xN-1)
        XeTVec(j+1) = 1i*(Lb - b*intCoord(j,0,xNb));
    end
end