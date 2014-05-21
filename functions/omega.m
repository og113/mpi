%function to find omega matrix of dl[7] - kuznetsov appendix
%output is omega matrix
%arguments are Nt and N
function Xomega = omega(xNt,xN)
    global L mass;
    xa = L/(xN-1);
    Xomega = complex(zeros(xN));
    h = zeros(xN);
    for j=1:xN
        if j>1 && j<xN
            h(j,j) = mass^2 + 2/xa^2;
            h(j,j+1) = -1/xa^2;
            h(j,j-1) = -1/xa^2;
        elseif j==1
            h(j,j) = mass^2 + 2/xa^2;
            h(j,j+1) = -1/xa^2;
        elseif j==xN
            h(j,j) = mass^2 + 2/xa^2;
            h(j,j-1) = -1/xa^2;
        end
    end
    
    [V,D] = eig(h);
    
    for j=1:xN
        for k=1:xN
            for l=1:xN
                Xomega(j,k) = Xomega(j,k) + xa*D(l,l)*V(j,l)*V(k,l);
            end
        end
    end
    
end