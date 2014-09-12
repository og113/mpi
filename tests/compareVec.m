%function to compare to files containing phi
%arguments are filenames 1 and 2
%outputs are vec1 phi2 and diff and vecname i.e. Cp or minusDS
function [vec1,vec2,diff] = compareVec(file1,file2, vecname)
    global Bdim Nb N;

    iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

    tempVec1 = loadVec(file1,vecname);
    tempVec2 = loadVec(file2,vecname);
    
    if (length(tempVec1)>Bdim)
        vec1 = vecComplex(tempVec1,Bdim);  
    else
        vec1 = tempVec1;
    end
    
    if (length(tempVec2)>Bdim)
        vec2 = vecComplex(tempVec2,Bdim);  
    else
        vec2 = tempVec2;
    end
    
    diff = vec1-vec2;
    
    diff(abs(diff)<1e-10) = 0; %removing small values
    
    t = imag(eTVec(Nb,N));
    x = xVec(Nb,N);
    
    subplot(3,2,1)
    plot3(t,x,real(vec1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(vec1)')
    subplot(3,2,2)
    plot3(t,x,imag(vec1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(vec1)')
    
    subplot(3,2,3)
    plot3(t,x,real(vec2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(vec2)')
    subplot(3,2,4)
    plot3(t,x,imag(vec2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(vec2)')
    
    subplot(3,2,5)
    plot3(t,x,real(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(diff)')
    subplot(3,2,6)
    plot3(t,x,imag(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(diff)')
    
end