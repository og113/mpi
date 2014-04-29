%function to give t(locNum+1)-t(locNum) in euclidean domain when this exists
%argument is locNum
%use is for links
function Xdt = dt(locNum)
    global Nt;
    if intCoord(locNum,0,Nt)==(Nt-1)
        Xdt=0;
    else
        Xdt = eCoord(locNum+1,0) - eCoord(locNum,0);
    end
end