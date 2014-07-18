%short function to plot phi output of script bubble
%arguments are p and N
function plotBubble(p,XN)
    XEdim = length(p)-1;
    p(XEdim+1) = [];
    XNt = XEdim/XN;
    t = rtVec(XEdim,XNt);
    x = xVec(XNt,XN);
    plot3(t,x,p,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('phi)')
end