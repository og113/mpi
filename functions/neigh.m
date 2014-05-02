%finds location number of neighbouring points in a lattice
%periodic in space but not time, gives -1 if no neighbour in given
%direction
%arguments are locNum, direction, sign and Nt
function Xneigh = neigh(locNum, direction, sign, xNt)
    global d N;
    Xneigh = -1; %//this is the result if there are no neighbours for the given values of the argument
    if direction==0
        if sign==1 && intCoord(locNum,0,xNt)~=(xNt-1)
            Xneigh = locNum+1;
        elseif sign==-1 && intCoord(locNum,0,xNt)~=0
            Xneigh = locNum-1;
        else
            disp('neigh error 1');
        end
    elseif intCoord(locNum,direction,xNt)==0 && sign==-1
        Xneigh = locNum + (N-1)*N^(direction-1)*xNt;
    elseif intCoord(locNum,direction,xNt)==(N-1) && sign==1
        Xneigh = locNum - (N-1)*N^(direction-1)*xNt;
    else
        Xneigh = locNum + sign*N^(direction-1)*xNt;
    end
end
    