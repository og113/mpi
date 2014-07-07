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
    subplot(1,2,1)
    plot3(t,x,real(tPhi),'x')
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('re(phi)')
    subplot(1,2,2)
    plot3(t,x,imag(tPhi),'x')
    xlabel('re(t)-im(t)'), ylabel('x'), zlabel('imag(phi)')
end