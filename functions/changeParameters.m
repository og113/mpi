%function to apply variation of parameters given user inputs
%arguments are: newParameter, parameterName and inputP
%where inputP is the input phi i.e. bubble, periodic instanton etc.
%NB - has not been checked for consistency.
function changeParameters(newParameter, parameterName, inputP)
    global d N Na Nb Nc NT L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon angle amp;
    if strcmp(parameterName,'mass')
        mass = newParameter;
        X = mass*R;
        epsilon = 2*mass^3/lambda/R/3;
        v =  mass/sqrt(lambda);
    elseif strcmp(parameterName,'epsilon')
        epsilon = newParameter;
        lambda = 2*mass^3/epsilon/R/3;
        v =  mass/sqrt(lambda);
    elseif strcmp(parameterName,'L')
        L = newParameter;
        a = L/(N-1);
    elseif strcmp(parameterName,'X')
        X = newParameter;
        mass = X/R;
        epsilon = 2*mass^3/lambda/R/3;
        v =  mass/sqrt(lambda);
    elseif strcmp(parameterName,'N')
        Na = floor(Na*(newParameter/N));
        Nb = floor(Nb*(newParameter/N));
        Nc = floor(Nc*(newParameter/N));
        N = newParameter;
        NT = N + Na + Nb + Nc;
        a = L/(N-1);
        b = Lb/(Nb-1);
        La = (Na-1)*b;
        Lc = (Nc-1)*b;
        Adim = N*Na;
        Bdim = N*Nb;
        Cdim = N*Nc;
        Tdim = Adim + Bdim + Cdim;
    elseif strcmp(parameterName,'lambda')
        epsilon = epsilon*(lambda/newParameter); %epsilon scales like 1/lambda
        lambda = newParameter;
        v =  mass/sqrt(lambda);
    elseif strcmp(parameterName,'amp')
        amp = newParameter; %amp is separate to everything else
    elseif strcmp(inputP,'b') || strcmp(inputP,'t') || strcmp(inputP,'f') %specific changes to be made for the spherical vacuum bubble or pure vacuum
        if strcmp(parameterName,'R')
            R = newParameter;
			L = 4*R;
			Lb = 2*R;
			X = mass*R;
			epsilon = 2*mass^3/lambda/R/3;
			v =  mass/sqrt(lambda);
			a = L/(N-1);
			b = Lb/(Nb-1);
            La = (Na-1)*b;
            Lc = (Nc-1)*b;
        elseif strcmp(parameterName,'Lb')
            Lb = newParameter;
            b = Lb/(Nb-1);
            La = (Na-1)*b;
            Lc = (Nc-1)*b;
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
            La = (Na-1)*b;
            X = mass*R;
			epsilon = 2*mass^3/lambda/R/3;
			v =  mass/sqrt(lambda);
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
            La = b*(Na-1);
            Lc = b*(Nc-1);
        end
    else
        disp('change parameters error');
    end