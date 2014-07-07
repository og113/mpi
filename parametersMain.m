%file generating derived global values of parameters for 'p' in main, in
%with input data struct from loaded file
function parameters(inputData)
    global d N Nt Ntm NT NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
    global R X lambda mass v epsilon angle theta;
    
    %main global parameters
    d = inputData.d;
    N = inputData.N;
    NtonN = inputData.NtonN;
    NtmonNt = inputData.NtmonNt;
    R = inputData.R;
    mass = inputData.mass;
    lambda = inputData.lambda;
    L = inputData.L;
    Lt = inputData.Lt;
    theta = 0;
    
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
    angle = asin(Lt/R/2);
    a = L/(N-1);
    b = Lt/(Nt-1);
    Ltm = (Ntm-1)*b;
end
        