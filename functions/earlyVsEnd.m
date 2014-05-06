%short function to plot early phi versus endphi
%arguments are early phi and end phi (as Cp)
function plotPhi(earlyPhi,endPhi)
    global Nt;
    subplot(1,2,1)
    x = xVec(Nt);
    t = imag(eTVec);
    plot3(t,x,real(earlyPhi),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(earlyPhi)');
    subplot(1,2,2)
    plot3(t,x,real(endPhi),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(endPhi)');
end