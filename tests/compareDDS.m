%function to compare DDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
%gives D1,D2,diff,maximum
function [XaDDS,XbDDS,XcompareDDS,Xmax] = compareDDS
    cutoff = eps;
    %mData = load('data/picEarly01.mat');
    
    %mData.DDSm(abs(mData.DDSv)<cutoff)=[]; %eliminating zeros
    %mData.DDSn(abs(mData.DDSv)<cutoff)=[];
    %mData.DDSv(abs(mData.DDSv)<cutoff)=[];
    
    %XmDDS = sparse(mData.DDSm,mData.DDSn,mData.DDSv);
    
    load data/D1.dat;
    D1(abs(D1(:,3))<cutoff,:)=[];
    XaDDS = spconvert(D1);
    
    load data/D2.dat;
    D2(abs(D2(:,3))<cutoff,:)=[];
    XbDDS = spconvert(D2);
    
    XcompareDDS = XaDDS - XbDDS;
    [i,j,s] = find(XcompareDDS);
    i(abs(s)<cutoff)=[];
    j(abs(s)<cutoff)=[];
    s(abs(s)<cutoff)=[];
    XcompareDDS = sparse(i,j,s);
    
    %subplot(1,3,1)
    %spy(XaDDS)
    %subplot(1,3,2)
    %spy(XbDDS)
    %subplot(1,3,3)
    %spy(XcompareDDS)
    
    tempV = nonzeros(XcompareDDS);
    Xmax = max(abs(tempV));
end

