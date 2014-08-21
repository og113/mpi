%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
function parameters(inputP,pot)
    global d N Na Nb Nc NT L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon dE minima angle theta;
    
    %main global parameters
    d = 2;
    N = 100;
    Na = 80; %changed later in 'p'
    Nb = 80;
    Nc = 32;
    lambda = 1/10;
    theta = 0;
    dE = 0.034;
%%%%%%%%%%%%%%%%%% - potentials
    clear x epsi;
    if pot==1
        v = 1.5811; %also a main global parameter
        epsilon = dE; %first guess
        eV = @(x,epsi) (lambda/8.0)*(x.^2-v^2)^2 - (epsi/2.0/v)*(x-v);
        edV = @(x,epsi) (x*lambda*(x.^2 - v^2))/2 - epsi/2.0/v;
        V = @(x) eV(x,epsilon);
        dV = @(x) edV(x,epsilon);
    elseif pot==2
        v = 0.4; %also a main global parameter
        epsilon = 0.75; %first guess
        W = @(x) exp(-x.^2)*(x + x.^3 + x.^5);
        dW = @(x) exp(-x^2)*(- 2*x^6 + 3*x^4 + x^2 + 1);
        eV = @(x,epsi) (1/2)*(x+1).^2*(1-epsi*W((x-1)/v));
        edV = @(x,epsi) (x+1)*(1-epsi*W((x-1)/v)) - (1/2)*(x+1).^2*(epsi/v)*dW((x-1)/v);
        V = @(x) eV(x,epsilon);
        dV = @(x) edV(x,epsilon);
    else
        disp('choice of potential not available, in parameters');
    end
%%%%%%%%%%%%%%%% - working out epsilon from dE
    test = 1;
    minCloseness = 1e-14;
    while test(end)>minCloseness
        oldMinima = fzero(dV,-3); %going from minus 3 to 3
        for j=1:50
            tempMinima = fzero(dV,-3+j*(6/50));
            if tempMinima>(oldMinima(end)+1e-15) || tempMinima<(oldMinima(end)-1e-15)
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
        for j=1:50
            tempMinima = fzero(dV,-3+j*(6/50));
            if tempMinima>(newMinima(end)+1e-15) || tempMinima<(newMinima(end)-1e-15)
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
        if length(test)>50
            disp('parameters looped over 50 times');
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
    mass =  v*sqrt(lambda);
    R = 2*mass^3/lambda/dE/3;
    X = mass*R;
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lb = 40;
		L = 3.8*R;
		a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
        if R>Lb || R>2*L
            disp('R is too large');
        end
    elseif inputP=='p' || inputP == 'q' || inputP == 'i'
        Lb = 40;
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