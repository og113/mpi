%function to give t(locNum+1)-t(locNum) in euclidean domain when this exists
%argument is locNum
%use is for links
function Xtdt = tdt(locNum)
    global NT;
    if intCoord(locNum,0,NT)==(NT-1)
        Xtdt=0;
    else
        Xtdt = coord(locNum+1,0) - coord(locNum-1,0);
    end
end