%function to give real coordinate value in euclidean domain
%arguments are locNUm and direction, ouput is scalar
function XreCoord = reCoord(locNum, direction, xNt)
    global Lb L a b;
    if direction==0
        XreCoord = Lb - b*intCoord(locNum,0,xNt);
    else
        XreCoord = -L/2 + a*intCoord(locNum,direction,xNt);
    end
end
