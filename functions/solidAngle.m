%function to calculate the solid angle in a given dimension
%argument is dimension
function XsolidAngle = solidAngle(dimension)
    if rem(dimension,2)==0
        XsolidAngle = 2*pi^(dimension/2)/factorial(dimension/2 - 1);
    elseif rem(dimension,2)==1
        XsolidAngle = 2^dimension*pi^((dimension-1)/2)*factorial((dimension-1)/2)/factorial(dimension-1);
    else
        disp('solid angle error');
    end
end