%function to test whether the newton-raphson loop has converged and to flag
%up various options if it hasn't
%ouput is [wait,aq]
%arguments are: runsCount, runsTest, aq, stopTime, action
function [Xwait, Xaq] = convergenceQuestions(aq, runsCount, stopTime, action, gaussianTest)
    global N Nb
    Xaq = aq;
    if runsCount > 1000
        disp('over 1000 runs without convergence - stopping n-r loop');
        Xwait = 0;
    elseif stopTime > 600
        disp(['time greater than ',num2str(stopTime),', number of n-r loops: ',num2str(runsCount) ]);
        Xwait = input('keep waiting?(1=y, 0=n) ');
        tic;
    elseif gaussianTest(end)>eps*2*N*Nb
        disp(['gaussian elimination didnt work well, norm of DDS*delta-minusDS = ',num2str(gaussianTest(end))]);
        Xwait = 0;
    else
        Xwait = 1;
    end
end