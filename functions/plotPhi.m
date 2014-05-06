%short function to plot real part of phi in euclidean domain
%argument is Cp
function plotPhi(Cp)
    global Nt;
    x = xVec(Nt);
    t = imag(eTVec);
    plot3(t,x,imag(Cp),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi)');