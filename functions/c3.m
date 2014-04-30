%like count except it counts mod 3 - an alternative name for count3
%that is the first three calls give 1, then the next three give 2 etc
%reset with 'clear c3'

function Xcount3 = count3
    persistent y;
    if isempty(y)
        y = 1;
    else
        y = y + 1;
    end
    Xcount3 = ceil(y/3);
end