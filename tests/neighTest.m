%function to test neigh, only works in 2d
%no arguments
function neighTest
    global Edim Nt
    for j=0:(Edim-1)
        vecToPrint{1} = ['(',num2str(intCoord(j,0,Nt)),',',num2str(intCoord(j,1,Nt)),'): '];
        for k=0:3
            sign = (-1)^k;
            direc = floor(k/2);
            if neigh(j,direc,sign,Nt)~=-1
                vecToPrint{k+2} = ['(',num2str(intCoord(neigh(j,direc,sign,Nt),0,Nt)),',',num2str(intCoord(neigh(j,direc,sign,Nt),1,Nt)),') '];
            else
                vecToPrint{k+2} = [];
            end
        end
        pause(0.1);
        for k=1:length(vecToPrint)
            fprintf('%12s',vecToPrint{k});
        end
            fprintf('\n');
    end
end