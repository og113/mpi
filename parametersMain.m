%file generating derived global values of parameters for 'p' in main, in
%with input data struct from loaded file
function parameters(inputData)
    global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim;
    global R epsilon dE theta;
    
    %main global parameters
    d = inputData.d;
    N = inputData.N;
    Na = inputData.Na;
    Nb = inputData.Nb;
    Nc = inputData.Nc;
    R = inputData.R;
    L = inputData.L;
    Lb = inputData.Lb;
    dE = inputData.dE;
    theta = 0;
    
    %derived global parameters
    NT = Na + Nb + Nc;
    Adim = N*Na;
    Bdim = N*Nb;
    Cdim = N*Nc;
    Tdim = N*NT;
    epsilon = inputData.epsilon;
    a = L/(N-1);
    b = Lb/(Nb-1);
    La = Na*b;
    Lc = Nc*b;
end
        