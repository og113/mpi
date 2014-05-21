%short function to plot phi output of script bubble
%argument is p
function plotBubble(p)
    XEdim = length(p)-1;
    XNt = sqrt(XEdim);
    t = rtVec(XEdim,XNt);
    x = xVec(XNt);
    if length(p)==(XEdim+1)
        p(XEdim+1) = [];
    end
    plot3(t,x,p,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('phi)')
end