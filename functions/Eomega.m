%function to find a variant of omega matrix of dl[7] - kuznetsov appendix
%output is variant of omega matrix - this is for evaluating E, there is an
%extra power of w_k over the usual omega matrix
%argument is N
function Xomega = omega(xN)
    global L;
    xa = L/(xN-1);
    Xomega = complex(zeros(xN));
    h = zeros(xN);
    for j=1:xN
        if j>1 && j<xN
            h(j,j) = 1 + 2/xa^2;
            h(j,j+1) = -1/xa^2;
            h(j,j-1) = -1/xa^2;
        elseif j==1
            h(j,j) = 1 + 2/xa^2;
            h(j,j+1) = -1/xa^2;
        elseif j==xN
            h(j,j) = 1 + 2/xa^2;
            h(j,j-1) = -1/xa^2;
        end
    end
    
    [V,D] = eig(h);
    
    for j=1:xN
        for k=1:xN
            for l=1:xN
                Xomega(j,k) = Xomega(j,k) + xa*V(j,l)*V(k,l)*D(l,l);
            end
        end
    end
    
end