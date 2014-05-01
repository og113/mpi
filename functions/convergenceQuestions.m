%function to test whether the newton-raphson loop has converged and to flag
%up various options if it hasn't
%ouput is [wait,aq]
%arguments are: runsCount, runsTest, aq, stopTime, action
function [Xwait,Xaq] = convergenceQuestions(runsCount, runsTest, aq, stopTime, action)
    if runsCount > 1000
        disp('over 1000 runs without convergence - stopping n-r loop');
        Xwait = 0;
    end
    if stopTime > 600
        disp(['time greater than ',num2str(stopTime),', number of n-r loops: ',num2str(runsCount) ]);
        printWait = input('print phi and action on next loop? (y/n)');
        if printWait == 'y'
            Xaq.printChoice = 'p';
            Xaq.printRun = runsCount + 1;
            disp(['action = ',num2str(action)]);
        end
        Xwait = input('keep waiting?(1=y, 0=n) ');
    end
end