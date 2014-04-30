%calculates coordinates in euclidean space
%argument is locNUm, ouput is vector
function XeCoords = eCoords(locNum)
    global d Nt Lt L a b;
    XeCoords = zeros(d,1);
    for l=1:d
        if l==1
            XeCoords(l) = 1i*(-Lt + b*intCoord(locNum,0,Nt));
        else
            XeCoords(l) = -L/2 + a*intCoord(locNum,l-1,Nt);
        end
    end
end