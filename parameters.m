%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
function parameters(inputP)
    global d N Na Nb Nc NT L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon angle theta;
    
    %main global parameters
    d = 2;
    N = 80;
    Na = 96;
    Nb = 80;
    Nc = 64;
    R = 32;
    mass = 1/2;
    lambda = 1/10;
    theta = 0;
    
    %derived global parameters
    NT = Na + Nb + Nc;
    Adim = Na*N;
    Bdim = Nb*N;
    Cdim = Nc*N;
    Tdim = NT*N;
    X = mass*R;
    epsilon = 2*mass^3/lambda/R/3;
    v =  mass/sqrt(lambda);
    
    %parameters specific to inputP
    if inputP=='b' || inputP=='f' || inputP=='t'
        Lb = 2*R;
		L = 3.5*R;
		a = L/(N-1);
		b = Lb/(Nb-1);
        La = (Na-1)*b;
        Lc = (Nc-1)*b;
    elseif inputP=='p'
        Lb = 1.0*R;
        angle = asin(Lb/R);
        Ltemp = 3*R;
        L = 3*Lb*tan(angle);
        a = L/(N-1);
		b = Lb/(Nb-1);
        La = (Na-1)*b;
        Lc = (Nc-1)*b;
    else
        disp('parameter error');
    end
end
        