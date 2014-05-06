%function to turn a 2*dimensional vector into a complex dimensional vector
%arguments are: the 2*dimensional vector, and the dimension
function XvecComplex = vecComplex(realVec, totDimension)
    XvecComplex = complex(zeros(totDimension,1));
    for j=0:(totDimension-1)
        XvecComplex(j+1) = realVec(2*j+1) + 1i*realVec(2*j+2);
    end
end