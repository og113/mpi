%file generating derived global values of parameters for 'p' in main, in
%with input data struct from loaded file
function parameters(inputData)
    global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim;
    global R X lambda mass v epsilon angle theta;
    
    %main global parameters
    d = inputData.d;
    N = inputData.N;
    Na = inputData.Na;
    Nb = inputData.Nb;
    Nc = inputData.Nc;
    R = inputData.R;
    mass = inputData.mass;
    lambda = inputData.lambda;
    L = inputData.L;
    Lb = inputData.Lb;
    theta = 0;
    
    %derived global parameters
    NT = Na + Nb + Nc;
    Adim = N*Na;
    Bdim = N*Nb;
    Cdim = N*Nc;
    Tdim = N*NT;
    X = mass*R;
    epsilon = 2*mass^3/lambda/R/3;
    v =  mass/sqrt(lambda);   
    angle = asin(Lb/R);
    a = L/(N-1);
    b = Lb/(Nb-1);
    La = Na*b;
    Lc = Nc*b;
end
        