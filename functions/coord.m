%function to give coordinate value on whole space
%arguments are locNUm and direction, ouput is scalar
function Xcoord = coord(locNum, direction)
    global NT Na Nb La Lb L a b;
    if direction==0
        t = intCoord(locNum,0,NT);
        if t<Na
            Xcoord =  -La+b*t + 1i*Lb;
        elseif t<(Na+Nb)
            Xcoord = 1i*(Lb - b*(t-Na));
        else
            Xcoord = (t-Na-Nb+1)*b;
        end
    else
        Xcoord = -L/2 + a*intCoord(locNum,direction,NT);
    end
end
