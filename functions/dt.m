%function to give t(locNum+1)-t(locNum) in euclidean domain when this exists
%argument is locNum
%use is for links
function Xbdt = bdt(locNum)
    global Nb;
    if intCoord(locNum,0,Nb)==(Nb-1)
        Xbdt=0;
    else
        Xbdt = eCoord(locNum+1,0) - eCoord(locNum-1,0);
    end
end