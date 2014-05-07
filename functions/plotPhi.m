%short function to plot real part of phi in euclidean domain
%argument is Cp
function plotPhi(Cp)
    global Nt;
    x = xVec(Nt);
    t = imag(eTVec);
    subplot(1,2,1)
    plot3(t,x,real(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi)')
end