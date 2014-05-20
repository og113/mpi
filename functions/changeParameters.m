%function to apply variation of parameters given user inputs
%arguments are: newParameter, parameterName and inputP
%where inputP is the input phi i.e. bubble, periodic instanton etc.
%NB - has not been checked for consistency.
function changeParameters(newParameter, parameterName, inputP)
    global d N Nt Ntm NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
    global R X lambda mass v epsilon theta;
    if parameterName == 'mass'
        mass = newParameter;
        X = mass*R;
        lambda = 2*(d-1)*mass^3/epsilon/R/3;
        v =  mass/sqrt(lambda);
    elseif parameterName == 'epsilon'
        epsilon = newParameter;
        lambda = 2*(d-1)*mass^3/epsilon/R/3;
        v =  mass/sqrt(lambda);
    elseif parameterName == 'L'
        L = newParameter;
        a = L/(N-1);
    elseif parameterName == 'X'
        X = newParameter;
        mass = X/R;
        lambda = 2*(d-1)*mass^3/epsilon/R/3;
        v =  mass/sqrt(lambda);
    elseif inputP == 'b' | inputP == 't' | inputP =='f' %specific changes to be made for the spherical vacuum bubble or pure vacuum
        Lt = 2*R;
		L = 4*R;
		a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
        if parameterName == 'N'
            N = newParameter;
			Nt = floor(N*NtonN);
            Ntm = floor(Nt*NtmonNt);
			a = L/(N-1);
			b = Lt/(Nt-1);
            Ltm = (Ntm-1)*b;
			Edim = N^(d-1)*Nt;
            Mdim = N^(d-1)*Ntm;
            Tdim = Edim + Mdim;
        elseif parameterName == 'R'
            R = newParameter;
			L = 4*R;
			Lt = 2*R;
			X = mass*R;
			lambda = 2*(d-1)*mass^3/epsilon/R/3.0;
			v =  mass/sqrt(lambda);
			a = L/(N-1);
			b = Lt/(Nt-1);
            Ltm = (Ntm-1)*b;
        elseif parameterName == 'Lt'
            Lt = newParameter;
            b = Lt/(Nt-1);
            Ltm = (Ntm-1)*b;
        elseif parameterName == 'lambda'
            epsilon = epsilon*(lambda/newParameter); %epsilon scales like one over lambda
            lambda = newParameter;
            R = 2*(d-1)*mass^3/epsilon/lambda/3;
            X = mass*R;
            v =  mass/sqrt(lambda);
            L = 4*R;
            Lt = 2*R;
            a = L/(N-1);
            b = Lt/(Nt-1);
            Ltm = b*(Ntm-1);
        end
    elseif inputP == 'p' %specific changes to be made for the periodic instanton
        Lt = 1.2*R/2; %Lt = T/2, where T is the period of the periodic instanton
        theta = asin(Lt/R);
        L = 3*Lt*tan(theta); %Lt*tan(theta) is thin-wall analytic width of periodic instanton bubble
        a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
        if parameterName == 'N'
            N = newParameter;
            Nt = floor(N*NtonN);
            Ntm = floor(Nt*NtmonNt);
            a = L/(N-1);
            b = Lt/(Nt-1);
            Ltm = b*(Ntm-1);
            Edim = N^(d-1)*Nt;
            Mdim = N^(d-1)*Ntm;
            Tdim = Edim + Mdim;
        elseif parameterName == 'R'
            R = newParameter;
            Lt = 1.2*R/2;
            theta = asin(Lt/R);
            L = 3*Lt*tan(theta);
            a = L/(N-1);
            b = Lt/(Nt-1);
            Ltm = (Ntm-1)*b;
            X = mass*R;
			lambda = 2*(d-1)*mass^3/epsilon/R/3.0;
			v =  mass/sqrt(lambda);
        elseif parameterName == 'Lt'
            Lt = newParameter;
            R = 2*Lt/1.2;
            theta = asin(Lt/R);
            L = 3*Lt*tan(theta);
            X = mass*R;
            lambda = 2*(d-1)*mass^3/epsilon/R/3.0;
			v =  mass/sqrt(lambda);
            a = L/(N-1);
            b = Lt/(Nt-1);
            Ltm = b*(Ntm-1);
        elseif parameterName == 'lambda'
            epsilon = epsilon*(lambda/newParameter); %epsilon scales like 1/lambda
            lambda = newParameter;
            R = 2*(d-1)*mass^3/epsilon/lambda/3;
            X = mass*R;
            v =  mass/sqrt(lambda);
            Lt = 1.2*R/2;
            theta = asin(Lt/R);
            L = 3*Lt*tan(theta);
            a = L/(N-1);
            b = Lt/(Nt-1);
            Ltm = b*(Ntm-1);
        end
    end