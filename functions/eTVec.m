%gives time in euclidean space as a vector
%arguments are Nt and N
function XeTVec = eTVec(xNt, xN)
    global Lt b;
    XeTVec = zeros(xNt*xN,1);
    for j=0:(xNt*xN-1)
        XeTVec(j+1) = 1i*(Lt - b*intCoord(j,0,xNt));
    end
end