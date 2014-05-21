%function to give real coordinate value in euclidean domain
%arguments are locNUm and direction, ouput is scalar
function XreCoord = reCoord(locNum, direction)
    global Nt Lt L a b;
    if direction==0
        XreCoord = Lt - b*intCoord(locNum,0,Nt);
    else
        XreCoord = -L/2 + a*intCoord(locNum,direction,Nt);
    end
end
