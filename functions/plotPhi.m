%short function to plot real part of phi in euclidean domain
%arguments are Cp Nt and N
function plotPhi(Cp,xNt,xN)
    x = xVec(xNt,xN);
    t = imag(eTVec(xNt,xN));
    subplot(1,2,1)
    plot3(t,x,real(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi)')
end