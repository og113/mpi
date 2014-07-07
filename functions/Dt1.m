%function to give, in euclidean domain, (t(locNum+1)-t(locNum-1))/2 in the bulk and the existent half of this at time boundaries
%argument is locNum
%use is for sites
function XDt = Dt(locNum)
global Nb;
    if intCoord(locNum,0,Nb)==0
        XDt = (eCoord(locNum+1,0) - eCoord(locNum,0))/2;
    elseif intCoord(locNum,0,Nb)==(Nb-1)
        XDt = (eCoord(locNum,0) - eCoord(locNum-1,0))/2;
    else
        XDt = (eCoord(locNum+1,0) - eCoord(locNum-1,0))/2;
    end
end