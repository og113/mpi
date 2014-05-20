%function to compare minusDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
function [XmminusDS,XcminusDS,XminusDS] = compareminusDS
    global Edim Nt;
    mData = load('data/minusDS.mat');
    XmminusDS = mData.minusDS;
    XmminusDS(Edim+1) = [];
    
    load data/earlyv.dat;
    XcminusDS = -earlyv(:,3);
    x = earlyv(:,2);
    t = earlyv(:,1);
    
    subplot(1,3,1)
    plot3(t,x,XmminusDS,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('mMinusDS')
    
    subplot(1,3,2)
    plot3(t,x,XcminusDS,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('cMinusDS')
    
    difference = XcminusDS-XmminusDS;
    
    subplot(1,3,3)
    plot3(t,x,difference,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('difference')

end

