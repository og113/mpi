%function to apply variation of parameters given user inputs
%arguments are: newParameter, parameterName and inputP
%where inputP is the input phi i.e. bubble, periodic instanton etc.
%NB - has not been checked for consistency.
function changeParameters(newParameter, parameterName, inputP)
    global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon angle;
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
    elseif strcmp(inputP,'b') || strcmp(inputP,'t') || strcmp(inputP,'f') %specific changes to be made for the spherical vacuum bubble or pure vacuum
        Lb = 2*R;
		L = 4*R;
		a = L/(N-1);
		b = Lb/(Nb-1);
        La = (Na-1)*b;
        Lc = (Nc-1)*b;
        if strcmp(parameterName,'N')
            NaonN = Na/N;
            NbonN = Nb/N;
            NconN = Nc/N;
            N = newParameter;
            Na = floor(N*NaonN);
            Nb = floor(N*NbonN);
			Nc = floor(N*NconN);
            NT = N + Na;
			a = L/(N-1);
			b = Lb/(Nb-1);
            La = (Na-1)*b;
            Lc = (Nc-1)*b;
            Adim = N*Na;
            Bdim = N*Nb;
            Cdim = N*Nc;
            Tdim = Adim + Bdim + Cdim;
        elseif strcmp(parameterName,'R')
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
        elseif strcmp(parameterName,'lambda')
            epsilon = epsilon*(lambda/newParameter); %epsilon scales like one over lambda
            lambda = newParameter;
            R = 2*mass^3/epsilon/lambda/3;
            X = mass*R;
            v =  mass/sqrt(lambda);
            L = 4*R;
            Lb = 2*R;
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = b*(Na-1);
            Lc = b*(Nc-1);
        end
    elseif strcmp(inputP,'p') %specific changes to be made for the periodic instanton
        Lb = 1.2*R/2; %Lt = T/2, where T is the period of the periodic instanton
        angle = asin(Lb/R);
        L = 3*Lb*tan(angle); %Lt*tan(angle) is thin-wall analytic width of periodic instanton bubble
        a = L/(N-1);
		b = Lb/(Nb-1);
        La = (Na-1)*b;
        Lc = (Nc-1)*b;
        if strcmp(parameterName,'N')
            NaonN = Na/N;
            NbonN = Nb/N;
            NconN = Nc/N;
            N = newParameter;
            Na = floor(N*NaonN);
            Nb = floor(N*NbonN);
			Nc = floor(N*NconN);
            NT = Nb + Na;
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = b*(Na-1);
            Lc = b*(Nc-1);
            Adim = N^(d-1)*Na;
            Bdim = N^(d-1)*Nb;
            Cdim = N^(d-1)*Nc;
            Tdim = Adim + Bdim + Cdim;
        elseif strcmp(parameterName,'R')
            R = newParameter;
            Lb = 1.2*R/2;
            angle = asin(Lb/R);
            L = 3*Lb*tan(angle);
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = (Na-1)*b;
            X = mass*R;
			epsilon = 2*mass^3/lambda/R/3;
			v =  mass/sqrt(lambda);
        elseif strcmp(parameterName,'Lb')
            Lb = newParameter;
            R = 2*Lb/1.2;
            angle = asin(Lb/R);
            L = 3*Lb*tan(angle);
            X = mass*R;
            epsilon = 2*mass^3/lambda/R/3;
			v =  mass/sqrt(lambda);
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = b*(Na-1);
            Lc = b*(Nc-1);
        elseif strcmp(parameterName,'lambda')
            epsilon = epsilon*(lambda/newParameter); %epsilon scales like 1/lambda
            lambda = newParameter;
            R = 2*mass^3/epsilon/lambda/3;
            X = mass*R;
            v =  mass/sqrt(lambda);
            Lb = 1.2*R/2;
            angle = asin(Lb/R);
            L = 3*Lb*tan(angle);
            a = L/(N-1);
            b = Lb/(Nb-1);
            La = b*(Na-1);
            Lc = b*(Nc-1);
        end
    end