%short function to plot phi output of script bubble
%argument is p
function plotBubble(p)
    global Nt Edim;
    t = rtVec;
    x = xVec(Nt);
    if length(p)==(Edim+1)
        p(Edim+1) = [];
    end
    plot3(t,x,p,'x')
    xlabel('im(t)'), ylabel('x'), zlabel('phi)')
end