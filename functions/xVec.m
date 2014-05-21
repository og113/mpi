%gives x as a vector
%argument is Nt or Ntm
function XxVec = xVec(xNt,xN)
    global L a;
    XxVec = zeros(xNt*xN,1);
    for j=0:(xNt*xN-1)
        XxVec(j+1) = -L/2 + a*intCoord(j,1,xNt);
    end
end