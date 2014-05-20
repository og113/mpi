%function to compare to files containing phi
%arguments are filenames 1 and 2
%outputs are vec1 phi2 and diff and vecname i.e. Cp or minusDS
function [vec1,vec2,diff] = comparePhi(file1,file2, vecname)
    global Edim Nt;

    iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

    vec1 = loadVec(file1,vecname);
    vec2 = loadVec(file2,vecname);
    
    checkLengths = @(y) iif( length(y)<Edim || length(y)>(Edim+1) ,@() error('lengths of vectors wrong'),true,@() 1);
    trimVector = @(y) iif( length(y)==(Edim+1) ,@() y([ones(1,Edim),0]),true, @() y);
    
    checkLengths(vec1);
    checkLengths(vec2);
    trimVector(vec1);
    trimVector(vec2);
    
    diff = vec1-vec2;
    
    diff(diff<1e-10) = 0; %removing small values
    
    t = imag(eTVec);
    x = xVec(Nt);
    
    subplot(3,2,1)
    plot3(t,x,real(vec1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(vec1)')
    subplot(3,2,2)
    plot3(t,x,imag(vec1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(vec1)')
    
    subplot(3,2,3)
    plot3(t,x,real(vec2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi2)')
    subplot(3,2,4)
    plot3(t,x,imag(vec2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi2)')
    
    subplot(3,2,5)
    plot3(t,x,real(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(diff)')
    subplot(3,2,6)
    plot3(t,x,imag(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(diff)')
    
end