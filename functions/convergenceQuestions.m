%function to test whether the newton-raphson loop has converged and to flag
%up various options if it hasn't
%ouput is [wait,aq]
%arguments are: runsCount, runsTest, aq, stopTime, action
function [Xwait, Xaq] = convergenceQuestions(aq, runsCount, stopTime, action)
    Xaq = aq;
    if runsCount > 1000
        disp('over 1000 runs without convergence - stopping n-r loop');
        Xwait = 0;
    elseif stopTime > 600
        disp(['time greater than ',num2str(stopTime),', number of n-r loops: ',num2str(runsCount) ]);
        printWait = input('print everything on following loops? (y/n) ','s');
        if printWait == 'y'
            Xaq.printChoice = 'e';
            Xaq.printRun = 0;
            disp(['action = ',num2str(action)]);
        end
        Xwait = input('keep waiting?(1=y, 0=n) ');
        tic;
    else
        Xwait = 1;
    end
end