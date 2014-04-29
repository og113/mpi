%calculates coordinates in euclidean space
%argument is locNUm, ouput is vector
function XeCoords = eCoords(locNum)
    global d Nt Lt L a b;
    XeCoords = zeros(d,1);
    for k=1:d
        if k==1
            XeCoords(k) = i*(-Lt + b*intCoord(locNum,0,Nt));
        else
            XeCoords(k) = -L/2 + a*intCoord(locNum,k-1,Nt);
        end
    end
end