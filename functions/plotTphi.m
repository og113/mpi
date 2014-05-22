%function to plot totPhi - may be worth subsuming this into plotPhi - by
%adding an argument to it
%arguments are tPhi, Nt, NT and N
function plotTphi(picOutStruct)
    global Lt;
    xNT = picOutStruct.N*picOutStruct.NtonN*(1 + picOutStruct.NtmonNt);
    xNt = picOutStruct.N*picOutStruct.NtonN;
    xN = picOutStruct.N;
    tPhi = picOutStruct.totalPhi;
    x = xVec(xNT,xN);
    t = real(tVec(xNt,xNT,xN)) + imag(tVec(xNt,xNT,xN)-1i*Lt);
    whos x
    whos t
    whos tPhi
    subplot(1,2,1)
    plot3(t,x,real(tPhi),'x')
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(tPhi),'x')
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('imag(phi)')
end