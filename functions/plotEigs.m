%quick function to look at negative eigenvectors of DDS
%arguments are V and lambda
function plotEigs(V,lambda)
    number = 12;
    for i=1:number
        subplot(2,number/2,i)
        plotBubble(V(:,i));
        title(['lambda = ',num2str(lambda(i))])
    end
end
    