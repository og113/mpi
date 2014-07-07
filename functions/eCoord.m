%calculates coordinates in euclidean space
%arguments are locNUm and direction, ouput is scalar
function XeCoord = eCoord(locNum, direction)
    global Nb Lb L a b;
    if direction==0
        XeCoord = 1i*(Lb - b*intCoord(locNum,0,Nb));
    else
        XeCoord = -L/2 + a*intCoord(locNum,direction,Nb);
    end
end