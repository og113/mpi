%function to compare to files containing phi
%arguments are filenames 1 and 2
%outputs are phi1 phi2 and diff
function [phi1,phi2,diff] = comparePhi(file1,file2)
    global Edim Nt;

    iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

    phi1 = loadPhi(file1);
    phi2 = loadPhi(file2);
    
    checkLengths = @(y) iif( length(y)<Edim || length(y)>(Edim+1) ,@() error('lengths of vectors wrong'),true,@() 1);
    trimVector = @(y) iif( length(y)==(Edim+1) ,@() y([ones(1,Edim),0]),true, @() y);
    
    checkLengths(phi1);
    checkLengths(phi2);
    trimVector(phi1);
    trimVector(phi2);
    
    diff = phi1-phi2;
    
    diff(diff<1e-10) = 0; %removing small values
    
    t = imag(eTVec);
    x = xVec(Nt);
    
    subplot(3,2,1)
    plot3(t,x,real(phi1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi1)')
    subplot(3,2,2)
    plot3(t,x,imag(phi1),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi1)')
    
    subplot(3,2,3)
    plot3(t,x,real(phi2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(phi2)')
    subplot(3,2,4)
    plot3(t,x,imag(phi2),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(phi2)')
    
    subplot(3,2,5)
    plot3(t,x,real(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('re(diff)')
    subplot(3,2,6)
    plot3(t,x,imag(diff),'x')
    xlabel('im(t)'), ylabel('x'), zlabel('imag(diff)')
    
end