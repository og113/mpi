%gives time in euclidean space as a vector
%no arguments
function XeTVec = eTVec
    global Nt Lt b Edim;
    XeTVec = zeros(Edim,1);
    for j=0:(Edim-1)
        XeTVec(j+1) = 1i*(Lt - b*intCoord(j,0,Nt));
    end
end