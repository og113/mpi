%function to compare DDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
function [XmDDS,XcDDS,XcompareDDS,Xmax] = compareDDS
    global Bdim;
    cutoff = 1e-13op    mData = load('data/DDS.mat');
    %mData.DDSm(mData.DDSn==(2*Bdim+1))=[]; %eliminating lagrange multiplier row and column
    %mData.DDSv(mData.DDSn==(2*Bdim+1))=[];
    %mData.DDSn(mData.DDSn==(2*Bdim+1))=[];
    %mData.DDSn(mData.DDSm==(2*Bdim+1))=[];
    %mData.DDSv(mData.DDSm==(2*Bdim+1))=[];
    %mData.DDSm(mData.DDSm==(2*Bdim+1))=[];
    
    mData.DDSm(mData.DDSv<cutoff)=[]; %eliminating zeros
    mData.DDSn(mData.DDSv<cutoff)=[];
    mData.DDSv(mData.DDSv<cutoff)=[];
    
    XmDDS = sparse(mData.DDSm,mData.DDSn,mData.DDSv);
    
    load data/cDDS.dat;
    cDDS(cDDS(:,3)<cutoff,:)=[];
    XcDDS = spconvert(cDDS);
    
    XcompareDDS = XmDDS - XcDDS;
    [i,j,s] = find(XcompareDDS);
    i(s<cutoff)=[];
    j(s<cutoff)=[];
    s(s<cutoff)=[];
    XcompareDDS = sparse(i,j,s);
    
    subplot(1,3,1)
    spy(XmDDS)
    subplot(1,3,2)
    spy(XcDDS)
    subplot(1,3,3)
    spy(XcompareDDS)
    
    tempV = nonzeros(XcompareDDS);
    Xmax = max(tempV);
end

