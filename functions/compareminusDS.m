%function to compare minusDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
function [XmminusDS,XcminusDS,XcompareminusDS] = compareminusDS
    global Edim;
    mData = load('data/minusDS.mat');
    XmminusDS = mData.minusDS;
    XmminusDS(2*Edim+1) = [];
    
    load data/cminusDS.dat;
    for j=0:(Edim-1)
        XcminusDS(2*j+1) = cminusDS(j+1,5);
        XcminusDS(2*j+2) = cminusDS(j+1,6);
    end
    XcminusDS = XcminusDS';
    
    XcompareminusDS = XmminusDS - XcminusDS;
    XcompareminusDS(XcompareminusDS<1e-10) = 0;
    
    subplot(1,3,1)
    plot(XmminusDS)
    subplot(1,3,2)
    plot(XcminusDS)
    subplot(1,3,3)
    plot(XcompareminusDS)
end

