%function to apply variation of parameters given user inputs
%arguments are: newParameter, parameterName and inputP
%where inputP is the input phi i.e. bubble, periodic instanton etc.
%NB - has not been checked for consistency.
%needs fixing post introduction of V1 and V2
function changeParameters(newParameter, parameterName, inputP, pot)
    global d N Na Nb Nc NT L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim;
    global R A epsilon dE theta angle amp minima;
    if strcmp(parameterName,'dE')
        dE = newParameter;
        [dE,epsilon,minima,S1] = potFn(pot,dE);
        R = S1/dE;
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
        if a>1 || b>1
            disp('lattice spacing too large');
        end
    elseif strcmp(parameterName,'L')
        L = newParameter;
        a = L/(N-1);
    elseif strcmp(parameterName,'N')
        Na = floor(Na*(newParameter/N));
        Nb = floor(Nb*(newParameter/N));
        Nc = floor(Nc*(newParameter/N));
        N = newParameter;
        NT = N + Na + Nb + Nc;
        a = L/(N-1);
        b = Lb/(Nb-1);
        La = Na*b;
        Lc = Nc*b;
        Adim = N*Na;
        Bdim = N*Nb;
        Cdim = N*Nc;
        Tdim = Adim + Bdim + Cdim;
    elseif strcmp(parameterName,'amp')
        amp = newParameter; %amp is separate to everything else
    elseif strcmp(inputP,'b') || strcmp(inputP,'t') || strcmp(inputP,'f') %specific changes to be made for the spherical vacuum bubble or pure vacuum
        if strcmp(parameterName,'R')
            R = newParameter;
			L = 4*R;
			Lb = 2*R;
			a = L/(N-1);
			b = Lb/(Nb-1);
            La = Na*b;
            Lc = Nc*b;
        elseif strcmp(parameterName,'Lb')
            Lb = newParameter;
            b = Lb/(Nb-1);
            La = Na*b;
            Lc = Nc*b;
        end
    elseif strcmp(inputP,'p') || strcmp(inputP,'q') || strcmp(inputP,'i') %specific changes to be made for the periodic instanton
        if strcmp(parameterName,'R')
            Lb = Lb*(newParameter/R);
            R = newParameter;
            angle = asin(Lb/R);
            L = 3*Lb*tan(angle);
            Ltemp = 4*R;
            if (L > abs(Ltemp) && Lb<=R) || (L < abs(Ltemp) && Lb>=R)%making sure to use the smaller of the two possible Ls
                L = Ltemp;
            end
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = Na*b;
            Lc = Nc*b;
        elseif strcmp(parameterName,'Lb')
            Lb = newParameter;
            angle = asin(Lb/R);
            L = 3*Lb*tan(angle);
            Ltemp = 4*R;
            if (L > Ltemp && Lb<=R) || (L < Ltemp && Lb>=R)%making sure to use the smaller of the two possible Ls
                L = Ltemp;
            end
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = b*Na;
            Lc = b*Nc;
        end
    else
        disp('change parameters error');
    end