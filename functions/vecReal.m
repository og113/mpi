%function to turn a 1*dimensional complex vector into a 2*dimensional real vector
%arguments are: the 1*dimensional vector, and the dimension e.g. Edim
function XvecReal = vecReal(cVec, totDimension)
    XvecReal = zeros(2*totDimension,1);
    for j=0:(totDimension-1)
        XvecReal(2*j+1) = real(cVec(j+1));
        XvecReal(2*j+2) = imag(cVec(j+1));
    end
end