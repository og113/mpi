%function to test whether the newton-raphson loop has converged and to flag
%- for main
%up various options if it hasn't
%ouput is [wait,aq]
%arguments are: runsCount, runsTest, stopTime, actionTest, vectorTest
function Xwait = convergenceQuestionsMain(runsCount, stopTime, actionTest, vectorTest)
    if runsCount > 1000
        disp('over 1000 runs without convergence - stopping n-r loop');
        Xwait = 0;
    elseif stopTime > 600
        disp(['time greater than ',num2str(stopTime),', number of n-r loops: ',num2str(runsCount) ]);
        Xwait = input('keep waiting?(1=y, 0=n) ');
        tic;
    elseif actionTest(end)>1e3 || vectorTest(end)>1e3
        disp(['wild nonconvergence, actionTest = ',num2str(actionTest(end)),', vectorTest = ',num2str(vectorTest(end))]);
        Xwait = 0;
    else
        Xwait = 1;
    end
end