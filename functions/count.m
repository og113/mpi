%function to count the number of calls of the function
%it takes no inputs and is persistent
%the first call is 1
%reset with 'clear count'
function Xcount = count
    persistent y;
    if isempty(y)
        y = 1;
    else
        y = y + 1;
    end
    Xcount = y;
end
    