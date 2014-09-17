%function to find epsilon etc given dE
%outputs are dE, epsilon, minima and S1
%arguments are pot and dE
function [OdE,Xepsilon,Xminima,S1] = potFn(pot,XdE)
    global A;
    clear x epsi;
    minCloseness = eps; %how close you have to get to XdE to end loop 
    if pot==1
        Xepsilon = XdE; %first guess
        eV = @(x,epsi) (1/8.0)*(x.^2-1.0).^2 - (epsi/2.0)*(x-1.0);
        edV = @(x,epsi) (x.*(x.^2 - 1.0))/2 - epsi/2.0;
        V = @(x) eV(x,Xepsilon);
        dV = @(x) edV(x,Xepsilon);
    elseif pot==2
        Xepsilon = 0.75; %first guess
        W = @(x) exp(-x.^2).*(x + x.^3 + x.^5);
        dW = @(x) exp(-x.^2).*(- 2*x.^6 + 3*x.^4 + x.^2 + 1);
        eV = @(x,epsi) (1/2)*(x+1).^2.*(1-epsi*W((x-1)/A));
        edV = @(x,epsi) (x+1).*(1-epsi*W((x-1)/A)) - (1/2)*(x+1).^2*(epsi/A)*dW((x-1)/A);
        V = @(x) eV(x,Xepsilon);
        dV = @(x) edV(x,Xepsilon);
    else
        disp('choice of potential not available, in parameters');
    end
%%%%%%%%%%%%%%%% - working out epsilon from dE
    test = 1;
    while test(end)>minCloseness
        oldMinima = fzero(dV,-3); %going from minus 3 to 3
        for j=1:100
            tempMinima = fzero(dV,-3+j*(6/100));
            addNew = 0;
            for k=1:length(oldMinima)
                if tempMinima>(oldMinima(k)+1e-10) || tempMinima<(oldMinima(k)-1e-12)
                    addNew = addNew + 1;
                end
            end
            if addNew==length(oldMinima)
                oldMinima(end+1) = tempMinima;
            end
        end
        if length(oldMinima)~=3
            disp('did not find all minima, try searching over a larger range, or check energy difference is small enough');
        end
        %olddE = abs(V(oldMinima(1))-V(oldMinima(3)));
        F = @(epsi) eV(oldMinima(1),epsi)-eV(oldMinima(3),epsi)-XdE;
        Xepsilon = fzero(F,Xepsilon);
        V = @(x) eV(x,Xepsilon);
        dV = @(x) edV(x,Xepsilon);
        newMinima = fzero(dV,-3); %going from minus 3 to 3
        for j=1:100
            tempMinima = fzero(dV,-3+j*(6/100));
            addNew = 0;
            for k=1:length(newMinima)
                if tempMinima>(newMinima(k)+1e-13) || tempMinima<(newMinima(k)-1e-13)
                    addNew = addNew + 1;
                end
            end
            if addNew==length(newMinima)
                newMinima(end+1) = tempMinima;
            end
        end
        if length(newMinima)~=3
            disp('did not find all minima, try searching over a larger range');
        end
        newdE = abs(V(newMinima(1))-V(newMinima(3)));
        OdE = newdE;
        for j=1:3
            testVec(j) = abs((newMinima(j)-oldMinima(j))/oldMinima(j));
        end
        testVec(4) = abs((newdE-XdE)/XdE);
        test(end+1) = max(testVec);
        if test(end)==test(end-1) && test(end)<1e-6
            break
        elseif length(test)>50
            disp('parameters looped over 50 times, consider increasing minCloseness');
        end
    end
    Xminima = newMinima;
    
    
    if pot == 1
        S1 = 2/3; %above isn't bad for pot==1 but 2/3 is a bit better
    else
        integrandS1 = @(x) (2.0*V(x)).^0.5;
        S1 = integral(integrandS1,Xminima(1)+3e-1,Xminima(3)-3e-1);
        if abs(imag(S1))>eps
            disp('S1 integral went too clse to boundaries');
        return
        end
    end
end
    