%calculates coordinates in euclidean space
%arguments are locNUm and direction, ouput is scalar
function XeCoord = eCoord(locNum, direction)
    global Nt Lt L a b;
    if direction==0
        XeCoord = 1i*(Lt - b*intCoord(locNum,0,Nt));
    else
        XeCoord = -L/2 + a*intCoord(locNum,direction,Nt);
    end
end