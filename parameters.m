%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
function parameters(inputP)
    global d N Nt Ntm NT NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
    global R X lambda mass v epsilon angle;
    
    %main global parameters
    d = 2;
    N = 64;
    NtonN = 1;
    NtmonNt = 2;
    R = 1e2;
    mass = 1;
    lambda = 1e-2;
    
    %derived global parameters
    Nt = NtonN*N;
    Ntm = NtmonNt*Nt;
    NT = Nt + Ntm;
    Edim = N^(d-1)*Nt;
    Mdim = N^(d-1)*Ntm;
    Tdim = Edim + Mdim;
    X = mass*R;
    epsilon = 2*(d-1)*mass^3/lambda/R/3;
    v =  mass/sqrt(lambda);
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lt = 2*R;
		L = 4*R;
		a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
    elseif inputP=='p'
        Lt = 1.2*R/2;
        angle = asin(Lt/R);
        L = 3*Lt*tan(angle);
        a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
    else
        disp('parameter error');
    end
end
        