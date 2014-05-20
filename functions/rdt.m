%function to give t(locNum+1)-t(locNum) in euclidean domain when this exists
%gives real dt, for program bubble
%argument is locNum
%use is for links
function Xrdt = rdt(locNum)
    global Nt;
    if intCoord(locNum,0,Nt)==(Nt-1)
        Xrdt=0;
    else
        Xrdt = reCoord(locNum+1,0) - reCoord(locNum,0);
    end
end