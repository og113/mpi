%function to compare DDS outputs from different programs/runs
%originally intended to compare c++ and matlab code
%no arguments
function [XmDDS,XcDDS,XcompareDDS] = compareDDS
    global Edim;
    mData = load('data/DDS.mat');
    mData.DDSm(mData.DDSn==(Edim+1))=[]; %eliminating lagrange multiplier row and column
    mData.DDSv(mData.DDSn==(Edim+1))=[];
    mData.DDSn(mData.DDSn==(Edim+1))=[];
    mData.DDSn(mData.DDSm==(Edim+1))=[];
    mData.DDSv(mData.DDSm==(Edim+1))=[];
    mData.DDSm(mData.DDSm==(Edim+1))=[];
    
    %mData.DDSm(mData.DDSv<1e-10)=[]; %eliminating zeros
    %mData.DDSn(mData.DDSv<1e-10)=[];
    %mData.DDSv(mData.DDSv<1e-10)=[];
    
    XmDDS = sparse(mData.DDSm,mData.DDSn,mData.DDSv);
    
    load data/cDDS.dat;
    %cDDS(cDDS(:,3)<1e-10,:)=[];
    XcDDS = spconvert(cDDS);
    
    XcompareDDS = XmDDS + XcDDS;
    [i,j,s] = find(XcompareDDS);
    i(abs(s)<1e-10)=[];
    j(abs(s)<1e-10)=[];
    s(abs(s)<1e-10)=[];
    XcompareDDS = sparse(i,j,s);
    
    subplot(1,3,1)
    spy(XmDDS)
    subplot(1,3,2)
    spy(XcDDS)
    subplot(1,3,3)
    spy(XcompareDDS)
end

