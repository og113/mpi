%gives x as a vector
%argument is Nt or Ntm
function XxVec = xVec(xNt)
    global L a Edim;
    XxVec = zeros(Edim,1);
    for j=0:(Edim-1)
        XxVec(j+1) = -L/2 + a*intCoord(j,1,xNt);
    end
end