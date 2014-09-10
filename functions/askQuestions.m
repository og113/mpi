%function to ask questions
%void function to change lots of arguments:
%inputP, perturbResponse, loopResponse, parameterChoice, minValue,
%maxValue, totalLoops, printChoice, printRun
function Xansw = askQuestions
    %Xansw.inputP = input('spherical bubble, periodic instanton, true vacuum or false vacuum? (b,p,t,f): ','s');
    %Xansw.perturbResponse = input('add small (10^-4) perturbations in real, imaginary or both parts of the field? (r/i/b/n): ','s');
    %Xansw.inputP = 'p';
    %Xansw.perturbResponse = 'n';
    %Xansw.pot = input('potential type 1 or 2? ');
    Xansw.loopResponse = input('loop through a parameter? (y/n): ','s');
    if Xansw.loopResponse == 'y'
        Xansw.parameterChoice = input('which parameter?: N, R, dE, L, Lb, amp: ' ,'s');
        Xansw.minValue = input('choose min value: ');
        Xansw.maxValue = input('choose max value: ');
        Xansw.totalLoops = input('choose number of loops: ');
    else
        Xansw.totalLoops = 1;
    end
end