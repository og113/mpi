%function to print action, phi, minusDS or DDS early
%arguments are: printChoice, printRun, runsCount and a variable argument giving the
%thing to print. for 'a' this should be a struct with the action, kinetic
%and potential terms in it
function earlyPrint(printChoice, printRun, runsCount, objectToPrint)
    if runsCount == printRun
        if printChoice == 'a' && isstruct(objectToPrint)
            disp(['kinetic = ',num2str(objectToPrint.kinetic)]);
            disp(['lambda potential = ',num2str(objectToPrint.potL)]);
            disp(['epsilon potential = ',num2str(objectToPrint.potE)]);
            disp(['action = ',num2str(objectToPrint.action)]);
        elseif printChoice=='p' && isvector(objectToPrint)
            save '/home/og113/Documents/mpi/data/phiEarly.dat' objectToPrint -ascii;
        elseif printChoice=='v' && isvector(objectToPrint)
            save '/home/og113/Documents/mpi/data/minusDS.dat' objectToPrint -ascii;
        elseif printChoice=='m' && ismatrix(objectToPrint)
            save '/home/og113/Documents/mpi/data/DDS.dat' objectToPrint -ascii;
        else
            disp('earlyPrint error');
        end
    end
end