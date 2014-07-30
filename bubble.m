%script to find the vacuum bubble 'b' using real phi everywhere, and
%dropping the complex parts
%the importance of this is to find the negative eigenvector and eigenvalue
%note that the zero mode constraint means that this doesn't work for 't' or
%'f'

%asking questions
inP = input('input bubble, true or false vacuum? (b,t,f) ','s');
loopResponse = input('loop through a parameter? (y/n): ','s');
if loopResponse == 'y'
    parameterChoice = input('which parameter?: N, R, mass, epsilon, L_t, X, lambda: ' ,'s');
    minValue = input('choose min value: ');
    maxValue = input('choose max value: ');
    totalLoops = input('choose number of loops: ');
else
    totalLoops = 1;
end
printChoice = input('print DDS matrix or -DS vector or the action (earlier) or phi or none? (m/v/a/p/n): ' ,'s');
if printChoice ~= 'n'
    printRun = input('choose run to print: ');
else
   printRun = -1; 
end

global d N NtonN L a b Bdim Nb Bdim Lb; %defining global variables
global R X lambda mass v epsilon;
parameters(inP); %assigning global variables according to parameters.m

tic; %starting the clock

for loop=0:(totalLoops-1) %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if loopResponse == 'y' && parameterChoice=='N'
        N = minValue + floor(loop*(maxValue - minValue)/(totalLoops-1));
        changeParameters(N,'N',inP);
    elseif loopResponse == 'y'
        loopParameter = minValue + loop*(maxValue - minValue)/(totalLoops-1);
        changeParameters (loopParameter,parameterChoice, inP);
    end
    
    S1 = 2*mass^3/3/lambda; %this is twice the value in the coleman paper
    twAction = -solidAngle(d)*epsilon*R^d/d + solidAngle(d)*R^(d-1)*S1; %thin-wall bubble action
    alpha = 8; %determines range over which tanh(x) is used
    action = 2;
    
    actionLast = 1; %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    runsTest = 1;
    closeness = 1e-4;
    minRuns = 6;
    
    p = zeros(Bdim+1,1); %phi, in the euclidean domain
    pNeg = zeros(2*Bdim,1); %negative eigenvector
    pZero = zeros(2*Bdim,1); %zero eigenvector
    
    syms y
    Roots = vpasolve(y^3 -v^2*y + epsilon/v/lambda,y); %solving V'(p)=0
    sort (Roots); %Roots sorted in ascending order
    
    for j=0:(Bdim-1) %assigning input phi according to inputP, note that the periodic instanton input is explicitly 2d
        rhoSqrd = 0;
        for k=0:(d-1)
            rhoSqrd = rhoSqrd + reCoord(j,k,Nb)^2;
        end
        rho = sqrt(rhoSqrd);
        if R<alpha/mass
            disp(['X = R*mass is too small. not possible to give thinwall input. it should be more that ',num2str(alpha)]);
        else
            if inP =='t'
                p(j+1) = Roots(1);
            elseif inP == 'f'
                p(j+1) = Roots(3);
            elseif inP == 'b'
                if rho<(R-alpha/mass)
                    p(j+1) = Roots(1);
                elseif rho>(R+alpha/mass)
                    p(j+1) = Roots(3);
                else
                    p(j+1) = v*tanh(mass*(rho-R)/2);
                    %%pNeg(2*j+1) = v/cosh(mass*(rho-R)/2)^2;
                end
            end
        end
    end
    
    p(Bdim+1) = v; %initializing Lagrange parameter for dp/dx zero mode
    
    %%if inP == 'b' || inP == 'p' %fixing norm of pNeg
        %%pNeg = pNeg/norm(pNeg);
    %%end
        
    Xwait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (runsTest > closeness || runsCount<minRuns) && Xwait%beginning newton-raphson loop
        runsTest = abs(action - actionLast)/abs(actionLast); %note this won't work if action goes to zero
        runsCount = runsCount + 1;
        actionLast = action;
        
        minusDS = zeros(Bdim+1,1); %-dS/d(p)
        %Chi0 = zeros(N,1); %to fix zero mode at final time, alla kuznetsov, dl[7], though not quite
        Chi1 = zeros(Bdim); %to fix zero mode everywhere

        %for j=1:(N-2) %explicitly 2d, evaluating Chi0, equal to the zero mode at t=(Nb-1)
        %    Chi0(j+1) = p((j+2)*Nb)-p(j*Nb); %shouldn't be forward or backward derivative as this will move boundary phi in a direction.
        %end
        %Chi0(N) = p(Nb)-p((N-1)*Nb);
        %Chi0(1) = p(2*Nb)-p(N*Nb);
        %if norm(Chi0)~=0
        %    Chi0 = Chi0/norm(Chi0);
        %end
        
        for j=0:(Bdim-1)
           Chi1(j+1) = p(neigh(j,1,1,Nb)+1) - p(neigh(j,1,-1,Nb)+1);
        end
        if norm(Chi1)~=0
            Chi1 = Chi1/norm(Chi1);
        end
        
        nonz = 5*(Bdim-2) + 5*N; %number of nonzero elements of DDS
        DDSm = zeros(nonz,1); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonz,1); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonz,1); %values of non-zero elements of DDS - don't forget to initialize DDS
        
        kinetic = 0; %initializing to zero
        potL = 0;
        potE = 0;
        clear c3;
        
        for j=0:(Bdim-1)
            t = intCoord(j,0,Nb);
            x = intCoord(j,1,Nb);

            if t==(Nb-1)
                siteMeasure = -a*b/2.0; %for sites in time
                linkMeasure = -a*b; %for links in time
            
                potL = potL + siteMeasure*(lambda/8)*(p(j+1)^2-v^2)^2;
                potE = potE + siteMeasure*epsilon*(p(j+1)-v)/v/2;
                if x~=(N-1) && x<(N-1)%evaluating spatial kinetic part
                    kinetic = kinetic + siteMeasure*(p(j+Nb+1)-p(j+1))^2/a^2/2;
                elseif x==(N-1) %avoinding using neigh and modulo as its slow
                    kinetic = kinetic + siteMeasure*(p(j+1-Nb*(N-1))-p(j+1))^2/a^2/2;
                end
            
                if inP=='b' %boundary conditions
                    DDSm(c3) = j+1; DDSn(c3) = j+1; DDSv(c3) = 1.0/b; %zero time derivative
                    DDSm(c3) = j+1; DDSn(c3) = j; DDSv(c3) = -1.0/b;
                elseif inP=='t' || inP=='f'
                    DDSm(c3) = j+1; DDSn(c3) = j+1; DDSv(c3) = 1.0; %zero change at boundary
                end
                
                %DDSm(c3) = j+1; DDSn(c3) = Bdim+1; DDSv(c3) = a*Chi0(x+1); %zero mode lagrange constraint
                %minusDS(j+1) = -a*p(Bdim+1)*Chi0(x+1);
                DDSm(c3) = j+1; DDSn(c3) = Bdim+1; DDSv(c3) = siteMeasure*Chi0(j+1); %zero mode lagrange constraint
                minusDS(j+1) = -siteMeasure*p(Bdim+1)*Chi1(j+1);
                minusDS(Bdim+1) = minusDS(Bdim+1) - siteMeasure*Chi1(j+1)*p((j+1)*Nb);        
                DDSm(c3) = Bdim+1; DDSn(c3) = (j+1)*Nb; DDSv(c3) = siteMeasure*Chi1(j+1); %at t=Nb-1
                
            else
                if t==0
                    siteMeasure = -a*b/2.0; %for sites in time
                    linkMeasure = -a*b; %for links in time
                    dtj = -b;

                    potL = potL + siteMeasure*(lambda/8.0)*(p(j+1)^2-v^2)^2;
                    potE = potE + siteMeasure*epsilon*(p(j+1)-v)/v/2;
                    if x~=(N-1) && x<(N-1)%evaluating spatial kinetic part
                        kinetic = kinetic + siteMeasure*(p(j+Nb+1)-p(j+1))^2/a^2/2;
                    elseif x==(N-1) %avoinding using neigh and modulo as its slow
                        kinetic = kinetic + siteMeasure*(p(j+1-Nb*(N-1))-p(j+1))^2/a^2/2;
                    end
                    kinetic = kinetic + linkMeasure*((p(j+2) - p(j+1))/dtj)^2/2;
                    
                    DDSm(c3) = j+1; DDSn(c3) = j+1; DDSv(c3) = 1; %zero change at boundary
                    
                    DDSm(c3) = j+1; DDSn(c3) = Bdim+1; DDSv(c3) = siteMeasure*Chi0(j+1); %zero mode lagrange constraint
                    minusDS(j+1) = -siteMeasure*p(Bdim+1)*Chi1(j+1);
                    minusDS(Bdim+1) = minusDS(Bdim+1) - siteMeasure*Chi1(j+1)*p((j+1)*Nb);        
                    DDSm(c3) = Bdim+1; DDSn(c3) = (j+1)*Nb; DDSv(c3) = siteMeasure*Chi1(j+1); %at t=Nb-1
                else
                    siteMeasure = -a*b; %for sites in time
                    linkMeasure = -a*b; %for links in time
                    dtj = -b;

                    potL = potL + siteMeasure*(lambda/8)*(p(j+1)^2-v^2)^2;
                    potE = potE + siteMeasure*epsilon*(p(j+1)-v)/v/2;
                    if x~=(N-1) && x<(N-1)%evaluating spatial kinetic part
                        kinetic = kinetic + siteMeasure*(p(j+Nb+1)-p(j+1))^2/a^2/2;
                    elseif x==(N-1) %avoinding using neigh and modulo as its slow
                        kinetic = kinetic + siteMeasure*(p(j+1-Nb*(N-1))-p(j+1))^2/a^2/2;
                    end
                    kinetic = kinetic + linkMeasure*((p(j+2) - p(j+1))/dtj)^2/2;
                    
                    for k=0:(2*d-1)
                        sign = (-1)^k;
                        %%deltaSign = (sign-1)/2; %deltaSign=0 if sign=+1 and deltaSign=-1 if sign=-1
                        direc = floor(k/2);
                        if direc == 0
                            minusDS(j+1) = minusDS(j+1) + a*p(j+sign+1)/dtj;
                            DDSm(c3) = j+1; DDSn(c3) = (j+sign)+1; DDSv(c3) = -a/dtj;
                        else
                            neighb = neigh(j,direc,sign,Nb);
                            minusDS(j+1) = minusDS(j+1) + siteMeasure*p(neighb+1)/a^2;
                            DDSm(c3) = j+1; DDSn(c3) = neighb+1; DDSv(c3) = -siteMeasure/a^2;
                        end
                    end                    
                    minusDS(j+1) = minusDS(j+1) - siteMeasure*(2.0*p(j+1)/a^2 + 2.0*p(j+1)/dtj^2 ...
                                    + (lambda/2.0)*p(j+1)*(p(j+1)^2-v^2) + epsilon/2.0/v);
                    DDSm(c3) = j+1; DDSn(c3) = j+1; DDSv(c3) = siteMeasure*( 2.0/a^2 + 2.0/dtj^2 + (lambda/2.0)*(3.0*p(j+1)^2-v^2) );
                    
                    DDSm(c3) = j+1; DDSn(c3) = Bdim+1; DDSv(c3) = siteMeasure*Chi0(j+1); %zero mode lagrange constraint
                    minusDS(j+1) = -siteMeasure*p(Bdim+1)*Chi1(j+1);
                    minusDS(Bdim+1) = minusDS(Bdim+1) - siteMeasure*Chi1(j+1)*p((j+1)*Nb);        
                    DDSm(c3) = Bdim+1; DDSn(c3) = j+1; DDSv(c3) = siteMeasure*Chi1(j+1); %at t=Nb-1
                end
            end
        end
        %for j=0:(Bdim-1) %adding last row, with lagrange multiplier terms
            %minusDS(Bdim+1) = minusDS(Bdim+1) - a*Chi0(j+1)*p((j+1)*Nb);        
            %DDSm(c3) = Bdim+1; DDSn(c3) = (j+1)*Nb; DDSv(c3) = a*Chi0(j+1); %at t=Nb-1        
        %end
        clear c3; %returning persistent output of c3 to 1
        %kinetic = 2*kinetic; potL = 2*potL; potE = 2*potE; %as we only calculated half of bubble in time direction
        action = 1i*(kinetic + potL + potE);
        DDSm(DDSv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
        DDSn(DDSv==0) = [];
        DDSv(DDSv==0) = [];
        DDS = sparse(DDSm,DDSn,DDSv,Bdim+1,Bdim+1);

        if printChoice~='n' && runsCount == printRun %printing early if asked for
            if printChoice == 'a'
                disp(['kinetic = ',num2str(kinetic)]);
                disp(['lambda potential = ',num2str(potL)]);
                disp(['epsilon potential = ',num2str(potE)]);
                disp(['action = ',num2str(action)]);
            elseif printChoice == 'p'
                save data/bubbleEarly.mat p;
                disp(['printed phi in data/bubbleEarly.mat on run ',num2str(runsCount)]);
            elseif printChoice == 'v'
                save data/minusDS.mat minusDS;
                disp(['printed minusDS in data/minusDS.mat on run ',num2str(runsCount)]);
            elseif printChoice == 'm'
                spy(DDS);
                pause(5);
                save data/DDS.mat DDSm DDSn DDSv;
                disp(['printed DDS in data/DDS.mat on run ',num2str(runsCount)]);
            else
                disp('early print error');
            end
        end
        
        small = norm(minusDS);
        smaller = small/(2*Bdim+1);
        cutoff = 1e-6;
        normed = 0;
        if smaller < cutoff
            disp(['small = ',num2str(small)]);
            break
        else
            normed = 1;
            minusDS = minusDS/small;
            DDS = DDS/small;
            disp(['small = ',num2str(small)]);
        end
        
        c = condest(DDS);
        fprintf('%12s','condest = ');
        fprintf('%12g\n',c);

        %[orderRow,orderCol,r,s,cc,rr] = dmperm(DDS); %preordering - gets vector order (and perhaps a second vector) - options are colamd, colperm and dmperm (which may produce 2 ordering vectors)
        orderRow = dmperm(DDS);
        
        setup.type = 'ilutp'; %preconditioning - incomplete LU factorization does not increase the number of non-zero elements in DDS - options are 'nofill', 'ilutp' and 'crout'
        setup.droptol = 1e-6; %drop tolerance is the minimum ratio of (off-diagonal) abs(U_ij) to norm(DDS(:j))
        setup.thresh = 0; %if 1 forces pivoting on diagonal, 0 to turn off
        setup.udiag = 0; %if 1 this replaces zeros in upper diagonal with droptol, 0 to turn off
        %[Lo,Up] = ilu(DDS,setup);%[Lo,Up] = ilu(DDS(orderRow,orderCol),setup);
        [Lo,Up] = ilu(DDS(orderRow,:),setup);

        tol = 1e-8; %tolerance for norm(Ax-b)/norm(b), consider increasing if procedure is slow
        maxit = 50; %max number of iterations
        delta = zeros(Bdim+1);
        %[delta(orderCol),flag,relres,iter,resvec] = lsqr(DDS(orderRow,orderCol),minusDS(orderRow),tol,maxit,Lo,Up); %finding solution iteratively. consider changing bicg to bicgstab, bicgstabl, cgs, gmres, lsqr, qmr or tfqmr 
        %[delta,flag,relres,iter,resvec] = lsqr(DDS,minusDS,tol,maxit,Lo,Up);
        [delta,flag,relres,iter,resvec] = lsqr(DDS(orderRow,:),minusDS(orderRow),tol,maxit,Lo,Up);
        if flag ~=0 %flag = 0 means bicg has converged, if getting non-zero flags, output and plot relres ([delta,flag] -> [delta,flag,relres])
            if flag == 1
                disp('linear solver interated maxit times but did not converge!');
                semilogy(1:length(resvec),resvec/norm(minusDS),'-o');
                xlabel('Iteration number');
                ylabel('Relative residual');
                k = input('find smallest k eigs: ');
                kEigs = eigs(DDS,k,1e-6); %1e-6 is the region around which eigenvalues will be found
                disp(kEigs);
                disp('press any key to break and resume');
                pause;
                break;
            elseif flag == 2
                disp('preconditioner M was ill-conditioned!');
            elseif flag == 3
                disp('linear solver stagnated (two consecutive iterates were the same)!');
            elseif flag == 4
                disp('one of the scalar quantities calculated during linear solver became too small or too large to continue computing!');
            end
        end
        fprintf('%12s','iter = ');
        fprintf('%12g\n',iter);
        fprintf('%12s','relres = ');
        fprintf('%12g\n',relres);
        
        if normed
            minusDS = minusDS*small;
            DDS = DDS*small;
        end

        
        if size(delta,1)==1
            p = p + delta'; %p -> p'
        else
            p = p + delta;
        end

        %pNeg and pZer plus log(det(DDS)) stuff
        stopTime = toc;
        
        if runsCount > 1000 %convergence questions
            disp('over 1000 runs without convergence - stopping n-r loop');
            wait = 0;
        elseif stopTime > 600
            disp(['time greater than ',num2str(stopTime),', number of n-r loops: ',num2str(runsCount) ]);
            printWait = input('print phi and action on next loop? (y/n) ','s');
            if printWait == 'y'
                printChoice = 'p';
                printRun = runsCount + 1;
                disp(['action = ',num2str(action)]);
            end
            wait = input('keep waiting?(1=y, 0=n) ');
            tic;
        else
            wait = 1;
        end
        
    end %closing newton-raphson loop

    %propagating solution back in minkowskian time
    %checking that converged
    
    if loop==0 %printing to terminal
        fprintf('%12s','time', 'runs','d','N','X','re(action)','im(action)'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('\n');
    end
    fprintf('%12g',toc,runsCount,d,N,X,real(action));
    fprintf('%12g\n',imag(action));
    
    actionOut = fopen('data/bubbleAction.dat','a'); %saving action etc to file
    fprintf(actionOut,'%12g',toc,runsCount,d,N,X,real(action));
    fprintf(actionOut,'%12g\n',imag(action));
    fclose(actionOut);
    
    save( ['data/bubble',num2str(loop),'.mat'], 'p', 'minusDS', 'DDS');%saving phi, DDS and minusDS to file
    
    tic;
    
end%closing parameter loop