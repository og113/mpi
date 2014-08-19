%file containing values of global parameters before effect of function
%changeParameters.m
%argument is inputP
function parameters(inputP)
    global d N Na Nb Nc NT L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon angle theta;
    
    %main global parameters
    d = 2;
    N = 100;
    Na = 80; %changed later in 'p'
    Nb = 80;
    Nc = 32;
    R = 26;
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
        Lb = 40;
		L = 3.8*R;
		a = L/(N-1);
		b = Lb/(Nb-1); %b section includes both corner points
        La = Na*b;
        Lc = Nc*b;
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