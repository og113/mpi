%function to ask questions
%void function to change lots of arguments:
%inputP, perturbResponse, loopResponse, parameterChoice, minValue,
%maxValue, totalLoops, printChoice, printRun
function Xansw = askQuestions
    Xansw.inputP = input('spherical bubble, periodic instanton, true vacuum or false vacuum? (b,p,t,f): ','s');
    Xansw.perturbResponse = input('add small (10^-4) perturbations in real, imaginary or both parts of the field? (r/i/b/n): ','s');
    Xansw.loopResponse = input('loop through a parameter? (y/n): ','s');
    if Xansw.loopResponse == 'y'
        Xansw.parameterChoice = input('which parameter?: N, R, mass, epsilon, L_t, X, lambda: ' ,'s');
        Xansw.minValue = input('choose min value: ');
        Xansw.maxValue = input('choose max value: ');
        Xansw.totalLoops = input('choose number of loops: ');
    end
    Xansw.printChoice = input('print DDS matrix or -DS vector or the action (earlier) or phi or none? (m/v/a/p/n): ' ,'s');
    if Xansw.printChoice ~= 'n'
        Xansw.printRun = input('choose run to print: ');
    end
end