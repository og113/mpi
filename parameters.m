%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
%parameters are now the barred parameters where powers of m have been used
%to make things dimensionless and an overall factor of 1/lambda has be
%factored out of the action
function parameters(inputP,pot)
    global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim;
    global R A epsilon dE minima angle theta;
    
    %main global parameters
    d = 2;
    N = 100;
    Na = 80; %changed later in 'p'
    Nb = 80;
    Nc = 32;
    theta = 0;
    dE = 0.05;
    A = 0.4; %only for pot2
%%%%%%%%%%%%%%%%%% - potentials
    clear x epsi;
    if pot==1
        epsilon = dE; %first guess
        eV = @(x,epsi) (1/8.0)*(x.^2-1.0)^2 - (epsi/2.0)*(x-1.0);
        edV = @(x,epsi) (x*(x.^2 - 1.0))/2 - epsi/2.0;
        V = @(x) eV(x,epsilon);
        dV = @(x) edV(x,epsilon);
    elseif pot==2
        epsilon = 0.75; %first guess
        W = @(x) exp(-x.^2).*(x + x.^3 + x.^5);
        dW = @(x) exp(-x.^2)*(- 2*x.^6 + 3*x.^4 + x.^2 + 1);
        eV = @(x,epsi) (1/2)*(x+1).^2*(1-epsi*W((x-1)/A));
        edV = @(x,epsi) (x+1).*(1-epsi*W((x-1)/A)) - (1/2)*(x+1).^2*(epsi/A)*dW((x-1)/A);
        V = @(x) eV(x,epsilon);
        dV = @(x) edV(x,epsilon);
    else
        disp('choice of potential not available, in parameters');
    end
%%%%%%%%%%%%%%%% - working out epsilon from dE
    test = 1;
    minCloseness = eps;
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
            disp('did not find all minima, try searching over a larger range');
        end
        %olddE = abs(V(oldMinima(1))-V(oldMinima(3)));
        F = @(epsi) eV(oldMinima(1),epsi)-eV(oldMinima(3),epsi)-dE;
        epsilon = fzero(F,epsilon);
        V = @(x) eV(x,epsilon);
        dV = @(x) edV(x,epsilon);
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
        for j=1:3
            testVec(j) = abs((newMinima(j)-oldMinima(j))/oldMinima(j));
        end
        testVec(4) = abs((newdE-dE)/dE);
        test(end+1) = max(testVec);
        if test(end)==test(end-1) && test(end)<1e-6
            dE = newdE;
            break
        elseif length(test)>50
            disp('parameters looped over 50 times, consider increasing minCloseness');
        end
    end
    minima = newMinima;
%%%%%%%%%%%%%%%%%%%%%    
    
    %simply derived global parameters
    NT = Na + Nb + Nc;
    Adim = Na*N;
    Bdim = Nb*N;
    Cdim = Nc*N;
    Tdim = NT*N;
    R = 2.0/dE/3.0;
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lb = 1.5*R;
		L = 3.5*R;
		a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
        if R>Lb || R>2*L
            disp('R is too large');
        end
    elseif inputP=='p' || inputP == 'q' || inputP == 'i'
        Lb = 1.5*R;
        angle = asin(Lb/R);
        Ltemp = 3.8*R;
        L = 1.5*(1.5*Lb*tan(angle)); %need L larger than La and Lc to fit close to light-like waves
        if (L > Ltemp && Lb<=R) || (L < Ltemp && Lb>=R)%making sure to use the smaller of the two possible Ls
            L = Ltemp;
        end
        a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
    else
        disp('parameter error');
    end
end