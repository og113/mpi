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
    Na = 64; %changed later in 'p'
    Nb = 64;
    Nc = 32;
    theta = 0;
    dE = 0.05;
    A = 0.4; %only for pot2    
    
    [dE,epsilon,minima,S1] = potFn(pot,dE);
    
    %simply derived global parameters
    NT = Na + Nb + Nc;
    Adim = Na*N;
    Bdim = Nb*N;
    Cdim = Nc*N;
    Tdim = NT*N;
    R = S1/dE;
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lb = 2.0*R;
		L = 4.0*R;
		a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
        if R>Lb || R>2*L
            disp('R is too large');
        end
    elseif inputP=='p' || inputP == 'q' || inputP == 'i'
        Lb = 1.3*R;
        L = 3.0*R;
        if Lb<R
            angle = asin(Lb/R);
            Ltemp = 3.0*(1.5*Lb*tan(angle)); %need L larger than La and Lc to fit close to light-like waves
            if (L > Ltemp && Lb<=R) || (L < Ltemp && Lb>=R)%making sure to use the smaller of the two possible Ls
                L = Ltemp;
            end
        else
            disp('Lb>=R');
        end
        a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
    else
        disp('parameter error');
    end
    if a>1 || b>1
        disp('lattice spacing too large');
    end
end