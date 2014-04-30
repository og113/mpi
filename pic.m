%script to find the periodic instanton of the stable phi^4 with negative
%mass

aq.inputP = []; %struct to hold answers to questions aq short for 'answers to questions'
aq.perturbResponse = [];
aq.loopResponse = [];
aq.parameterChoice = [];
aq.minValue = [];
aq.maxvalue = [];
aq.totalLoops = [];
aq.printChoice = [];
aq.printRun = [];

aq = askQuestions; %asking questions

global d N Nt Ntm NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim; %defining global variables
global R X lambda mass v epsilon theta;
parameters(aq.inputP); %assigning global variables according to parameters.m

tic; %starting the clock

for loop=1:aq.totalLoops %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if aq.loopResponse == 'y' && aq.parameterChoice=='N'
        N = aq.minValue + floor(loop*(aq.maxValue - aq.minValue),(aq.totalLoops-1));
        changeParameters(N,'N',aq.inputP);
    elseif aq.loopResponse == 'y'
        loopParameter = aq.minValue + loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1);
        changeParameters (loopParameter,aq.parameterChoice, aq.inputP);
    end
    
    S1 = 2*mass^3/3/lambda; %this is twice the value in the coleman paper
    twAction = -solidAngle(d)*epsilon*R^d/d + solidAngle(d)*R^(d-1)*S1; %thin-wall bubble action
    alpha = 20; %determines range over which tanh(x) is used
    action = complex(0);
    
    actionLast = complex(1); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    runsTest = 1;
    closeness = 1e-4;
    minRuns = 6;
    
    p_e = zeros(2*Edim,1); %phi, in the euclidean domain
    perturbReal = zeros(Edim,1); %perturbations
    perturbImag = zeros(Edim,1);
    %pNeg = zeros(2*Edim,1); %negative eigenvector
    
    syms x
    roots = solve(x^3 -v^2*x + epsilon/v/lambda == 0); %solving V'(p)=0
    sort (roots); %roots sorted in ascending order

    if aq.perturbResponse ~='n' %assigning values to perturbations if user wants perturbations
        perturbReal = v*1e-4*rand(Edim,1);
        perturbImag = v*1e-4*rand(Edim,1);
        for j=1:(d-1)
            for k=0:(N-1)
                if inputP == 't' || inputP =='f' %pure vacua
                    perturbReal(k*Nt*N^(j-1)+1) = 0; %zero perturbation at both (euclidean) time boundaries
                    perturbReal((k+1)*Nt*N^(j-1)) = 0;
                    perturbImag(k*Nt*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nt*N^(j-1)) = 0;
                elseif inputP == 'b' %bubble
                    perturbReal(k*Nt*N^(j-1)+1) = 0; %zero perturbation at initial time
                    perturbReal((k+1)*Nt*N^(j-1)) = perturbReal((k+1)*Nt*N^(j-1)-1); %zero derivative at final time
                    perturbImag(k*Nt*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nt*N^(j-1)) = perturbImag((k+1)*Nt*N^(j-1)-1);
                elseif inputP == 'p' %periodic instanton
                    perturbReal(k*Nt*N^(j-1)+1) = perturbReal(k*Nt*N^(j-1)+2); %zero time derivative at both time boundaries
                    perturbReal((k+1)*Nt*N^(j-1)) = perturbReal((k+1)*Nt*N^(j-1)-1);
                    perturbImag(k*Nt*N^(j-1)+1) = perturbImag(k*Nt*N^(j-1)+2);
                    perturbImag((k+1)*Nt*N^(j-1)) = perturbImag((k+1)*Nt*N^(j-1)-1);
                end
            end
        end
    end
    
    for j=0:(Edim-1) %assigning input phi according to inputP, note that the periodic instanton input is explicitly 2d
        rhoSqrd = -eCoord(j,0)^2;
        rho1Sqrd = -(eCoord(j,0)+1i*Lt/2)^2; %explicitly 2d here as would need rho3 and rho4 for the 3d case etc.
        rho2Sqrd = -(eCoord(j,0)+1i*Lt/2)^2; 
        for k=1:(d-1)
            rhoSqrd = rhoSqrd + coord(j,k)^2;
            rho1Sqrd = rho1Sqrd + (coord(j,k)+R*cos(theta))^2;
            rho2Sqrd = rho2Sqrd + (coord(j,k)-R*cos(theta))^2;
        end
        rho = real(sqrt(rhoSqrd)); %rho should be real even without the real()
        rho1 = real(sqrt(rho1Sqrd));
        rho2 = real(sqrt(rho2Sqrd));
        if R<alpha/mass
            disp(['X = R*mass is too small. not possible to give thinwall input. it should be less that ',num2str(alpha)]);
        else
            if inputP =='t'
                p_e(2*j+1) = roots(1);
            elseif inputP == 'f'
                p_e(2*j+1) = roots(3);
            elseif inputP == 'b'
                if rho<(R-alpha/mass)
                    p_e(2*j+1) = roots(1);
                elseif rho>(R+alpha/mass)
                    p_e(2*j+1) = roots(3);
                else
                    p_e(2*j+1) = v*tanh(mass*(rho-R)/2);
                    %pNeg(2*j+1) = v/cosh(mass*(rho-R)/2)^2;
                end
            elseif inputP == 'p'
                if rho1<(R-alpha/mass) && rho2<(R-alpha/mass)
                    p_e(2*j+1) = roots(1);
                elseif rho1>(R+alpha/mass) || rho2>(R+alpha/mass)
                    p_e(2*j+1) = roots(3);
                elseif real(eCoord(j,1))<0 %explicitly 2d here, note that the coord should be real
                    p_e(2*j+1) = v*tanh(mass*(rho1-R)/2);
                    %pNeg(2*j+1) = v/cosh(mass*(rho1-R)/2)^2;
                elseif real(eCoord(j,1))>0
                    p_e(2*j+1) = v*tanh(mass*(rho2-R)/2);
                    %pNeg(2*j+1) = v/cosh(mass*(rho2-R)/2)^2;
                else
                    p_e(2*j+1) = roots(3); %if eCoord(j,1) == 0
                end
            end
            if aq.perturbResponse == 'r' || aq.perturbResponse =='b'
                p_e(2*j+1) = p_e(2*j+1) + perturbReal(j+1);
            end
            if aq.perturbResponse == 'i' || aq.perturbResponse =='b'
                p_e(2*j+2) = p_e(2*j+2) + perturbImag(j+1);
            end
        end
    end
    
    if inputP == 'p' %fixing input periodic instanton to have zero time derivative at time boundaries
        for j=0:(N-1)
            p_e(2*(j*Nt+1)+1) = (p_e(2*j*Nt+1) + p_e(2*(j*Nt+1)+1))/2;
            p_e(2*j*Nt+1) = (p_e(2*j*Nt+1) + p_e(2*(j*Nt+1)+1))/2;
            p_e(2*(j*Nt+1)+2) = (p_e(2*j*Nt+2) + p_e(2*(j*Nt+1)+2))/2;
            p_e(2*j*Nt+2) = (p_e(2*j*Nt+2) + p_e(2*(j*Nt+1)+2))/2;
            p_e(2*((j+1)*Nt-1)) = (p_e(2*((j+1)*Nt)) + p_e(2*((j+1)*Nt-1)))/2;
            p_e(2*((j+1)*Nt)) = (p_e(2*((j+1)*Nt)) + p_e(2*((j+1)*Nt-1)))/2;
            p_e(2*((j+1)*Nt-2)+2) = (p_e(2*((j+1)*Nt-1)+2) + p_e(2*((j+1)*Nt-2)+2))/2;
            p_e(2*((j+1)*Nt-1)+2) = (p_e(2*((j+1)*Nt-1)+2) + p_e(2*((j+1)*Nt-2)+2))/2;
        end
    end
    
    %if inputP == 'b' | inputP == 'p' %fixing norm of pNeg
        %norm = dot(pNeg,pNeg);
        %pNeg = pNeg/norm;
    %end
        
    wait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (runsTest > closeness || runsCount<minRuns) && wait %beginning newton-raphson loop
        runsTest = abs(action - actionLast)/abs(actionLast); %note this won't work if action goes to zero
        runsCount = runsCount + 1;
        actionLast = action;
        
        minusDS = zeros(2*Edim,1); %-dS/d(p_e)
        Cp_e = complex(zeros(Edim,1)); %complex p_e
        delta = zeros(2*Edim,1); %p'-p for newton-raphson loop
        pZero = zeros(2*Edim,1); %zero mode = d(p_e)/dx - note that in more than 2d we will need more pZeros
        nonZ = 0; %number of nonzero elements of DDS
        if inputP == 'b' || inputP == 't' || inputP == 'f'
            nonZ = 5 + 4*(Nt-2)*(2*d+1);
        elseif inputP == 'p'
            nonZ = 6 + 4*(Nt-2)*(2*d+1);
        else %number of nonzero elements when fixing zero mode and with full boundary conditions
            nonZ = 3*N^(d-1)+1 + 4*(Nt-2)*(2*d+1);
        end
        DDSm = zeros(nonZ); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonZ); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonZ); %values of non-zero elements of DDS
        
        action = complex(0); %initializing to zero
        kinetic = complex(0);
        potL = complex(0);
        potE = complex(0);
        
        for j=0:Edim
            t = intCoord(j,0,Nt);
            siteMeasure = a^(d-1)*Dt(j); %for sites in time
            linkMeasure = a^(d-1)*dt(j); %for links in time
            
            potL = potL - siteMeasure*(lambda/8)*((p_e(2*j+1)+1i*p_e(2*j+2))^2-v^2)^2;
            potE = potE - siteMeasure*epsilon*(p_e(2*j+1)+1i*p_e(2*j+2)-v)/v/2;
            for k=1:(d-1) %evaluating spatial kinetic part
                if intCoord(j,k,Nt)~=(N-1)
                    kinetic = kinetic - siteMeasure*(p_e(2*neigh.......
        
