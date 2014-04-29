%script to find the periodic instanton of the stable phi^4 with negative
%mass

%defining global variables
global d N Nt;
d = 2; %dimension
N = 64; %number of points in spatial directions
Nt = 64; %number of points in time direction

%struct to hold answers to questions
answ.inputP;
answ.perturbResponse;
answ.loopResponse;
answ.parameterChoice;
answ.minValue;
answ.maxvalue;
answ.totalLoops;
answ.printChoice;
answ.printRun;

answ = askQuestions;