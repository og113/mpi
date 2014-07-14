%script to find the periodic instanton of the stable phi^4 with negative
%mass

aq.inputP = ['b']; %struct to hold answers to questions aq short for 'answers to questions' - defaults in initialization
aq.perturbResponse = ['n'];
aq.loopResponse = ['y'];
aq.parameterChoice = ['N'];
aq.minValue = [32];
aq.maxValue = [64];
aq.totalLoops = [5];
aq.printChoice = ['n'];
aq.printRun = [];

%%%aq = askQuestions; %asking questions
inP = aq.inputP; %just because I write this a lot

global d N Nt Ntm NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim; %defining global variables
global R X lambda mass v epsilon theta;
parameters(inP); %assigning global variables according to parameters.m

tic; %starting the clock

for loop=1:aq.totalLoops %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if aq.loopResponse == 'y' && aq.parameterChoice=='N'
        N = aq.minValue + floor(loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1));
        changeParameters(N,'N',inP);
    elseif aq.loopResponse == 'y'
        loopParameter = aq.minValue + loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1);
        changeParameters (loopParameter,aq.parameterChoice, inP);
    end
    
    S1 = 2*mass^3/3/lambda; %this is twice the value in the coleman paper
    twAction = -solidAngle(d)*epsilon*R^d/d + solidAngle(d)*R^(d-1)*S1; %thin-wall bubble action
    alpha = 15; %determines range over which tanh(x) is used
    action = complex(2);
    
    actionLast = complex(1); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    runsTest = 1;
    closeness = 1e-4;
    minRuns = 3;
    
    p = zeros(2*Edim,1); %phi, in the euclidean domain
    perturbReal = zeros(Edim,1); %perturbations
    perturbImag = zeros(Edim,1);
    %%pNeg = zeros(2*Edim,1); %negative eigenvector
    
    syms x
    roots = vpasolve(x^3 -v^2*x + epsilon/v/lambda,x); %solving V'(p)=0
    sort (roots); %roots sorted in ascending order

    if aq.perturbResponse ~='n' %assigning values to perturbations if user wants perturbations
        perturbReal = v*1e-4*rand(Edim,1);
        perturbImag = v*1e-4*rand(Edim,1);
        for j=1:(d-1)
            for k=0:(N-1)
                if inP == 't' || inP =='f' %pure vacua
                    perturbReal(k*Nt*N^(j-1)+1) = 0; %zero perturbation at both (euclidean) time boundaries
                    perturbReal((k+1)*Nt*N^(j-1)) = 0;
                    perturbImag(k*Nt*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nt*N^(j-1)) = 0;
                elseif inP == 'b' %bubble
                    perturbReal(k*Nt*N^(j-1)+1) = 0; %zero perturbation at initial time
                    perturbReal((k+1)*Nt*N^(j-1)) = perturbReal((k+1)*Nt*N^(j-1)-1); %zero derivative at final time
                    perturbImag(k*Nt*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nt*N^(j-1)) = 0;
                elseif inP == 'p' %periodic instanton
                    perturbReal(k*Nt*N^(j-1)+1) = perturbReal(k*Nt*N^(j-1)+2); %zero time derivative at both time boundaries
                    perturbReal((k+1)*Nt*N^(j-1)) = perturbReal((k+1)*Nt*N^(j-1)-1);
                    perturbImag(k*Nt*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nt*N^(j-1)) = 0;
                end
            end
        end
    end
    
    for j=0:(Edim-1) %assigning input phi according to inputP, note that the periodic instanton input is explicitly 2d
        t = eCoord(j,0);
        x = eCoord(j,1); %explicitly 2d
        rhoSqrd = -t^2 + x^2;
        rho1Sqrd = -(t-1i*Lt/2)^2 + (x+R*cos(theta))^2; %explicitly 2d here as would need rho3 and rho4 for the 3d case etc.
        rho2Sqrd = -(t-1i*Lt/2)^2 + (x-R*cos(theta))^2; 
        rho = real(sqrt(rhoSqrd)); %rho should be real even without the real()
        rho1 = real(sqrt(rho1Sqrd));
        rho2 = real(sqrt(rho2Sqrd));
        if R<alpha/mass
            disp(['X = R*mass is too small. not possible to give thinwall input. it should be less that ',num2str(alpha)]);
        else
            if inP =='t'
                p(2*j+1) = roots(1);
            elseif inP == 'f'
                p(2*j+1) = roots(3);
            elseif inP == 'b'
                if rho<(R-alpha/mass)
                    p(2*j+1) = roots(1);
                elseif rho>(R+alpha/mass)
                    p(2*j+1) = roots(3);
                else
                    p(2*j+1) = v*tanh(mass*(rho-R)/2);
                    %%pNeg(2*j+1) = v/cosh(mass*(rho-R)/2)^2;
                end
            elseif inP == 'p'
                if rho1<(R-alpha/mass) && rho2<(R-alpha/mass)
                    p(2*j+1) = roots(1);
                elseif rho1>(R+alpha/mass) || rho2>(R+alpha/mass)
                    p(2*j+1) = roots(3);
                elseif real(x)>0 %explicitly 2d here, note that the coord should be real
                    p(2*j+1) = v*tanh(mass*(rho1-R)/2);
                    %%pNeg(2*j+1) = v/cosh(mass*(rho1-R)/2)^2;
                elseif real(x)<0
                    p(2*j+1) = v*tanh(mass*(rho2-R)/2);
                    %%pNeg(2*j+1) = v/cosh(mass*(rho2-R)/2)^2;
                else
                    p(2*j+1) = roots(1); %if eCoord(j,1) == 0
                end
            end
            if aq.perturbResponse == 'r' || aq.perturbResponse =='b'
                p(2*j+1) = p(2*j+1) + perturbReal(j+1);
            end
            if aq.perturbResponse == 'i' || aq.perturbResponse =='b'
                p(2*j+2) = p(2*j+2) + perturbImag(j+1);
            end
        end
    end
    
    if inP == 'p' %fixing input periodic instanton to have zero time derivative at time boundaries
        for j=0:(N-1)
            p(2*(j*Nt+1)+1) = (p(2*j*Nt+1) + p(2*(j*Nt+1)+1))/2;
            p(2*j*Nt+1) = (p(2*j*Nt+1) + p(2*(j*Nt+1)+1))/2;
            p(2*(j*Nt+1)+2) = (p(2*j*Nt+2) + p(2*(j*Nt+1)+2))/2;
            p(2*j*Nt+2) = (p(2*j*Nt+2) + p(2*(j*Nt+1)+2))/2;
            p(2*((j+1)*Nt-1)) = (p(2*((j+1)*Nt)) + p(2*((j+1)*Nt-1)))/2;
            p(2*((j+1)*Nt)) = (p(2*((j+1)*Nt)) + p(2*((j+1)*Nt-1)))/2;
            p(2*((j+1)*Nt-2)+2) = (p(2*((j+1)*Nt-1)+2) + p(2*((j+1)*Nt-2)+2))/2;
            p(2*((j+1)*Nt-1)+2) = (p(2*((j+1)*Nt-1)+2) + p(2*((j+1)*Nt-2)+2))/2;
        end
    end
    
    %%if inP == 'b' | inP == 'p' %fixing norm of pNeg
        %%pNeg = pNeg/norm(pNeg);
    %%end
        
    Xwait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (runsTest > closeness || runsCount<minRuns) && Xwait%beginning newton-raphson loop
        runsTest = abs(action - actionLast)/abs(actionLast); %note this won't work if action goes to zero
        runsCount = runsCount + 1;
        actionLast = action;
        
        minusDS = zeros(2*Edim,1); %-dS/d(p)
        Cp = vecComplex(p,Edim);
        %pZero = zeros(2*Edim,1); %zero mode = d(p)/dx - note that in more than 2d we will need more pZeros
        
        nonz = 0; %number of nonzero elements of DDS
        if inP == 'b' || inP == 't' || inP == 'f'
            nonz = 5*N^(d-1) + 4*N^(d-1)*(Nt-2)*(2*d+1);
        elseif inP == 'p'
            nonz = 6*N^(d-1) + 4*N^(d-1)*(Nt-2)*(2*d+1);
        else %number of nonzero elements when fixing zero mode and with full boundary conditions
            nonz = 5*Edim*N^(d-1) + N^(d-1) + 4*(Nt-2)*N^(d-1)*(2*d+1);
        end
        DDSm = zeros(nonz,1); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonz,1); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonz,1); %values of non-zero elements of DDS - don't forget to initialize DDS
        
        action = complex(0); %initializing to zero
        kinetic = complex(0);
        potL = complex(0);
        potE = complex(0);
        clear c3;
        
        for j=0:(Edim-1)
            t = intCoord(j,0,Nt);
            dtj = dt(j);
            dtjm = dt(j-1);
            Dtj = Dt(j);
            siteMeasure = a^(d-1)*Dtj; %for sites in time
            
            potL = potL - siteMeasure*(lambda/8)*(Cp(j+1)^2-v^2)^2;
            potE = potE - siteMeasure*epsilon*(Cp(j+1)-v)/v/2;
            for k=1:(d-1) %evaluating spatial kinetic part
                if intCoord(j,k,Nt)~=(N-1)
                    kinetic = kinetic - siteMeasure*(Cp(neigh(j,k,1,Nt)+1)-Cp(j+1))^2/a^2/2;
                end
            end
            if t==(Nt-1)
                if inP=='p' || inP=='b'
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1/b; %zero time derivative %don't forget to clear c3
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j-1)+1; DDSv(c3) = -1/b;
                elseif inP=='t' || inP=='f'
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1; %zero change at boundary
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1;
                end
            else
                kinetic = kinetic + a*(Cp(j+2) - Cp(j+1))^2/dtj/2;
                if t==0
                    if inP=='b' || inP=='t' || inP=='f'
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1; %zero change at boundary
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1;
                    elseif inP=='p'
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = -1/b; %zero time derivative
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+1)+1; DDSv(c3) = 1/b;
                    end
                else
                    for k=0:(2*d-1)
                        sign = (-1)^k;
                        dtd = dtj;
                        if sign==-1
                           dtd = dtjm; 
                        end
                        direc = floor(k/2);
                        if direc == 0
                            minusDS(2*j+1) = minusDS(2*j+1) + real(a*Cp(j+sign+1)/dtd);
                            minusDS(2*j+2) = minusDS(2*j+2) + imag(a*Cp(j+sign+1)/dtd);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+sign)+1; DDSv(c3) = -real(a/dtd);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+sign)+2; DDSv(c3) = imag(a/dtd);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*(j+sign)+1; DDSv(c3) = -imag(a/dtd);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*(j+sign)+2; DDSv(c3) = -real(a/dtd);
                        else
                            neighb = neigh(j,direc,sign,Nt);
                            minusDS(2*j+1) = minusDS(2*j+1) - real(Dtj*Cp(neighb+1)/a); %note Dt(j) = Dt(j+Nt) etc.
                            minusDS(2*j+2) = minusDS(2*j+2) - imag(Dtj*Cp(neighb+1)/a);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*+1; DDSv(c3) = real(Dtj/a);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*neighb+2; DDSv(c3) = -imag(Dtj/a);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+1; DDSv(c3) = imag(Dtj/a);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+2; DDSv(c3) = real(Dtj/a);
                        end
                    end
                    temp0 = a^(d-1)*(1/dtj + 1/dtjm);
                    temp1 = siteMeasure*(2*(d-1)*Cp(j+1)/a^2 + (lambda/2)*Cp(j+1)*(Cp(j+1)^2-v^2) + epsilon/2/v);
                    temp2 = siteMeasure*(2*(d-1)/a^2 + (lambda/2)*(3*Cp(j+1)^2 - v^2));

                    minusDS(2*j+1) = minusDS(2*j+1) + real(temp1 - temp0*Cp(j+1));
                    minusDS(2*j+2) = minusDS(2*j+2) + imag(temp1 - temp0*Cp(j+1));
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = real(-temp2 + temp0);
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+2; DDSv(c3) = imag(temp2 - temp0);
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+1; DDSv(c3) = imag(-temp2 + temp0);
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = real(-temp2 + temp0);
                end
            end
        end
        clear c3; %returning persistent output of c3 to 1
        kinetic = 2*kinetic; potL = 2*potL; potE = 2*potE; %as we only calculated half of bubble in time direction
        action = kinetic + potL + potE;
        DDSm(DDSv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
        DDSn(DDSv==0) = [];
        DDSv(DDSv==0) = [];
        DDS = sparse(DDSm,DDSn,DDSv,2*Edim,2*Edim);

        if aq.printChoice~='n' && runsCount == aq.printRun %printing early if asked for
            if aq.printChoice == 'a'
                disp(['kinetic = ',num2str(kinetic)]);
                disp(['lambda potential = ',num2str(potL)]);
                disp(['epsilon potential = ',num2str(potE)]);
                disp(['action = ',num2str(action)]);
            elseif aq.printChoice == 'p'
                t = eTVec;
                x = xVec(Nt);
                save data/phiEarly.mat t x Cp;
                disp(['printed phi in data/phiEarly.mat on run ',num2str(runsCount)]);
            elseif aq.printChoice == 'v'
                save data/minusDS.mat minusDS;
                disp(['printed minusDS in data/minusDS.mat on run ',num2str(runsCount)]);
            elseif aq.printChoice == 'm'
                spy(DDS);
                save data/DDS.mat DDS;
                disp(['printed DDS in data/DDS.mat on run ',num2str(runsCount)]);
            else
                disp('early print error');
            end
        end

        setup.type = 'ilutp'; %incomplete LU factorization does not increase the number of non-zero elements in DDS
        setup.droptol = 1e-6; %drop tolerance is the minimum ratio of (off-diagonal) abs(U_ij) to norm(DDS(:j))
        [Lo,Up] = ilu(DDS,setup);

        tol = 1e-12; %tolerance for norm(Ax-b)/norm(b), consider increasing if procedure is slow
        maxit = 20; %max number of iterations
        [delta,flag,relres,iter,resvec] = cgs(DDS,minusDS,tol,maxit,Lo,Up); %finding solution iteratively. consider changing bicg to bicgstab, bicgstabl, cgs, gmres, lsqr, qmr or tfqmr 
        if flag ~=0 %flag = 0 means bicg has converged, if getting non-zero flags, output and plot relres ([delta,flag] -> [delta,flag,relres])
            if flag == 1
                disp('cgs interated maxit times but did not converge!');
                semilogy(0:maxit,resvec/norm(minusDS),'-o');
                xlabel('Iteration number');
                ylabel('Relative residual');
                k = input('find smallest k eigs: ');
                kEigs = eigs(DDS,k,0);
                disp(kEigs);
                disp('press any key to break and resume');
                pause;
                break;
            elseif flag == 2
                disp('preconditioner M was ill-conditioned!');
            elseif flag == 3
                disp('cgs stagnated (two consecutive iterates were the same)');
            elseif flag == 4
                disp('one of the scalar quantities calculated during cgs became too small or too large to continue computing');
            end
        end

        p = p + delta; %p -> p'

        %pNeg and pZer plus log(det(DDS)) stuff

        stopTime = toc;
        [Xwait,aq] = convergenceQuestions(runsCount, runsTest, aq, stopTime, action); %discovering whether or not n-r has converged, and stopping if it is wildly out

    end %closing newton-raphson loop

    %propagating solution back in minkowskian time
    %checking that converged
    
    if loop==1 %printing to terminal
        fprintf('%12s','time', 'runs','d','N','X','re(action)','im(action)'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('\n');
    end
    fprintf('%12g',toc,runsCount,d,N,X,real(action));
    fprintf('%12g\n',imag(action));
    
    actionOut = fopen('data/picAction.dat','a'); %saving action to file
    fprintf(actionOut,'%12g',toc,runsCount,d,N,X,real(action));
    fprintf(actionOut,'%12g\n',imag(action));
    fclose(actionOut);
    
    t = eTVec; %saving phi to file
    x = xVec(Nt);
    phiOut = ['data/phi',num2str(loop),'.mat'];
    save(phiOut,'t','x','Cp');
    save data/phiEarly.mat t x Cp;
    
end%closing parameterChange loop
