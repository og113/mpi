%function to give, in euclidean domain, (t(locNum+1)-t(locNum-1))/2 in the bulk and the existent half of this at time boundaries
%gives rDt for function bubble
%argument is locNum
%use is for sites
function XrDt = rDt(locNum)
global Nt;
    if intCoord(locNum,0,Nt)==0
        XrDt = (reCoord(locNum+1,0) - reCoord(locNum,0))/2;
    elseif intCoord(locNum,0,Nt)==(Nt-1)
        XrDt = (reCoord(locNum,0) - reCoord(locNum-1,0))/2;
    else
        XrDt = (reCoord(locNum+1,0) - reCoord(locNum-1,0))/2;
    end
end