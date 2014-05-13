%function to plot totPhi - may be worth subsuming this into plotPhi - by
%adding an argument to it
function plotTphi(tPhi)
    global Nt Ntm Lt;
    x = xVec(Nt+Ntm);
    t = real(tVec) - imag(tVec-1i*Lt);
    subplot(1,2,1)
    plot3(t,x,real(tPhi),'x')
    xlabel('re(t)-im(t-1i*Lt)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(tPhi),'x')
    xlabel('re(t)-im(t-1i*Lt)'), ylabel('x'), zlabel('imag(phi)')
end