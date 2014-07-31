%short function to plot real part of phi in euclidean domain
%arguments are Cp Nt and N, can also put p instead of Cp
function plotPhi(Cp,xNt,xN)
    x = real(xVec(xNt,xN));
    t = imag(eTVec(xNt,xN));
    if length(Cp)>xNt*xN
        Cp = vecComplex(Cp,xNt*xN);
    end
    subplot(1,2,1)
    plot3(t,x,real(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi)')
    
    
end