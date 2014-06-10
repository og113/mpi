%function to give coordinate value on whole space
%arguments are locNUm and direction, ouput is scalar
function Xcoord = coord(locNum, direction)
    global NT Ntm Lt L a b;
    if direction==0
        t = intCoord(locNum,0,NT);
        if t<Ntm
            Xcoord =  b*(t-Ntm+1) + 1i*Lt;
        else
            Xcoord = 1i*(Lt - b*(t-Ntm));
        end
    else
        Xcoord = -L/2 + a*intCoord(locNum,direction,NT);
    end
end
