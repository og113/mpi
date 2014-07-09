%function to plot totPhi - may be worth subsuming this into plotPhi - by
%adding an argument to it
%argument is picOutStruct
function plotTphi(picOutStruct)
    global Lb;
    xNa = picOutStruct.Na;
    xNb = picOutStruct.Nb;
    xNc = picOutStruct.Nc;
    xNT = xNa + xNb + xNc;
    xN = picOutStruct.N;
    tPhi = picOutStruct.tCp;
    x = xVec(xNT,xN);
    t = real(tVec(xN,xNa,xNb,xNc)) + imag(-tVec(xN,xNa,xNb,xNc)+1i*Lb);
    eLogical = (t>0) & (t<Lb);
    ePhi = tPhi(eLogical);
    ex = x(eLogical);
    et = t(eLogical);
    mPhi = tPhi(~eLogical);
    mx = x(~eLogical);
    mt = t(~eLogical);
    subplot(1,2,1)
    plot3(et,ex,real(ePhi),'rx');
    hold;
    plot3(mt,mx,real(mPhi),'x');
    hold;
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(et,ex,imag(ePhi),'rx');
    hold;
    plot3(mt,mx,imag(mPhi),'x');
    hold;
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('imag(phi)')
end