%function to test dt
function dtTest
    global Edim Nt;
    for j=0:(Edim-1)
        vecToPrint{1} = ['(',num2str(intCoord(j,0,Nt)),',',num2str(intCoord(j,1,Nt)),'): '];
        vecToPrint{2} = [num2str(real(Dt(j))),' + ',num2str(imag(Dt(j))),'*i '];
        pause(0.1);
        for k=1:length(vecToPrint)
            fprintf('%12s',vecToPrint{k});
        end
            fprintf('\n');
    end
end