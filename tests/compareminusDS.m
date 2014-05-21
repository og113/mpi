%function to compare minusDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
function [XmminusDS,XcminusDS,XminusDS] = compareminusDS
    global Edim Nt;
    mData = load('data/minusDS.mat');
    XmminusDS = mData.minusDS;
    XmminusDS(2*Edim+1) = [];
    Cm = vecComplex(XmminusDS,Edim);
    
    load data/cminusDS.dat;
    XcminusDS = cminusDS(:,5)' + 1i*cminusDS(:,6)';
    
    load data/minusDS.dat;
    XminusDS = minusDS(:,3)';
    
    x = xVec(Nt);
    t = imag(eTVec);
    
    subplot(2,3,1)
    plot3(t,x,real(XcminusDS),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(XcminusDS)')
    subplot(2,3,4)
    plot3(t,x,imag(XcminusDS),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(XcminusDS)')
    
    subplot(2,3,2)
    plot3(t,x,real(Cm),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(Cm)')
    subplot(2,3,5)
    plot3(t,x,imag(Cm),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(Cm)')
    
    subplot(2,3,3)
    plot3(t,x,real(XminusDS),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(C)')
    subplot(2,3,6)
    plot3(t,x,imag(XminusDS),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(C)')

end

