%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
function parameters(inputP)
    global d N Nt Ntm NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
    global R X lambda mass v epsilon theta;
    
    %main global parameters
    d = 2;
    N = 64;
    NtonN = 1;
    NtmonNt = 1;
    R = 100;
    mass = 1;
    epsilon = 1e-4;
    
    %derived global parameters
    Nt = NtonN*N;
    Ntm = NtmonNt*Nt;
    Edim = N^(d-1)*Nt;
    Mdim = N^(d-1)*Ntm;
    Tdim = Edim + Mdim;
    X = mass*R;
    lambda = 2*(d-1)*mass^3/epsilon/R/3;
    v =  mass/sqrt(lambda);
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lt = 2*R;
		L = 4*R;
		a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
    elseif inputP=='p'
        Lt = 1.2*R;
        theta = asin(Lt/2/R);
        L = 1.5*Lt*tan(theta);
        a = L/(N-1);
		b = Lt/(Nt-1);
        Ltm = (Ntm-1)*b;
    else
        disp('parameter error');
    end
end
        