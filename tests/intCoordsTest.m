%script to test intCoords.m
%should clear workspace before use

paren = @(x, varargin) x(varargin{:});
%defining an anonymous function to take the components of a vector valued
%function (see
%http://blogs.mathworks.com/loren/2013/01/24/introduction-to-functional-programming-with-anonymous-functions-part-2/?s_eid=PSM_3365)
%for more information)

%defining global variables
global d N Nt;
d = 2; %dimension
N = 4; %number of points in spatial directions
Nt = 4; %number of points in time direction

%printing global variables
toDisplay = [ 'd = ', num2str(d), ', N = ', num2str(N), ', Nt = ', num2str(Nt)];
disp(toDisplay);

%asking form of test
displayAns = input('display single coordinate or whole array (s,w) ', 's');

%testing
if displayAns == 's'
    locNum = input ( 'give location: ');
    coords = zeros(d,1);
    coords = intCoords(locNum, Nt);
    disp(coords);
elseif displayAns == 'w'
    locNumVec = 0:(N^(d-1)*Nt-1);
    coordVec = zeros(N^(d-1)*Nt,1);
    for j=1:d
        for i= 1:(N^(d-1)*Nt-1)
            coordVec(i+1) = paren(intCoords(i,Nt),j);
        end
        subplot(d,1,j)
        plotTitle = ['direction = ',num2str(j)];
        plot(locNumVec,coordVec), title(plotTitle)
    end
end