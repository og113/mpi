%finds location number of neighbouring points in a lattice
%periodic in space but not time, gives -1 if no neighbour in given
%direction
%arguments are locNum, direction, sign and Nt
function Xneigh = neigh(locNum, direction, sign, Nt)
    global d N;
    Xneigh = -1; %//this is the result if there are no neighbours for the given values of the argument
    if direction==0
        if sign==1 && intCoord(locNum,0,Nt)~=(Nt-1)
            Xneigh = locNum+1;
        elseif sign==-1 && intCoord(locNum,0,Nt)~=0
            Xneigh = locNum-1;
        else
            disp('neigh error 1');
        end
    elseif intCoord(locNum,direction,Nt)==0 && sign==-1
        Xneigh = locNum + (N-1)*N^(direction-1)*Nt;
    elseif intCoord(locNum,direction,Nt)==(N-1) && sign==1
        Xneigh = locNum - (N-1)*N^(direction-1)*Nt;
    else
        Xneigh = locNum + sign*N^(direction-1)*Nt;
    end
end
    