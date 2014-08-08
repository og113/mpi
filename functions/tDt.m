%function to give, in euclidean domain, (t(locNum+1)-t(locNum-1))/2 in the bulk and the existent half of this at time boundaries
%argument is locNum
%use is for sites
function XtDt = tDt(locNum)
global NT;
    Xt = intCoord(locNum,0,NT);
    if Xt==0
        XtDt = (coord(locNum+1,0) - coord(locNum,0))/2;
    elseif Xt==(NT-1)
        XtDt = (coord(locNum,0) - coord(locNum-1,0))/2;
    else
        XtDt = (coord(locNum+1,0) - coord(locNum-1,0))/2;
    end
end