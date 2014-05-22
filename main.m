%script to solve boundary value problem with non-zero T and angle. input is
%the periodic instanton and then steps increasing angle or T.

global d N Nt Ntm NT NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim; %defining global variables
global R X lambda mass v epsilon angle theta;

fileNo = input('which data/picOut#.mat file to load? (#) '); %loading periodic instanton
data = load(['data/picOut',num2str(fileNo),'.mat']);
data.DDS = []; data.minusDS = []; data.Cp = []; %freeing some memory

disp(['Lt = T/2 = ', num2str(data.Lt)); %asking input questions
disp('angle = 0');
maxTheta = input('input final value of angle (input 0 for no steps) ');
%%maxLt = input('and input final value of Lt');
totalLoops = 1;
if maxTheta~=0
    totalLoops = input('and number of steps to get there ');
end

parametersMain(data); %assigning global variables according to data and parametersMain.m

p = data.totalPhi; data.totalPhi = []; %assigning phi from pic.m output, and freeing up memory

tic; %starting the clock

for loop=0:(totalLoops-1) %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if totalLoops>1
        theta = theta + loop*maxTheta/(totalLoops - 1);
    end
    
    alpha = 20; %determines range over which tanh(x) is used
    action = complex(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    actionLast = complex(1); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    runsTest = 1;
    closeness = 1e-4;
    minRuns = 6;
    
    %%pNeg = zeros(2*Tdim,1); %negative eigenvector
    %%pZero = zeros(2*Tdim,1); %zero eigenvector
    
    syms x
    roots = vpasolve(x^3 -v^2*x + epsilon/v/lambda,x); %solving V'(p)=0
    sort (roots); %roots sorted in ascending order

    if loop==0
        for j=0:(N-1) %fixing input periodic instanton to obey boundary conditions
            p(2*(j*NT+1)+1) = (p(2*j*NT+1) + p(2*(j*NT+1)+1))/2; %d(re(p))/dt=0 at initial time
            p(2*j*NT+1) = (p(2*j*NT+1) + p(2*(j*NT+1)+1))/2; %d(re(p))/dt=0 at initial time
            p(2*j*NT+2) = 0; %im(p)=0 at initial time
            p(2*((j+1)*NT-1)+1) = (p(2*((j+1)*NT-1)+1) + p(2*((j+1)*NT-2)+1))/2; %d(re(p))/dt=0 at final time
            p(2*((j+1)*NT-2)+1) = (p(2*((j+1)*NT-1)+1) + p(2*((j+1)*NT-2)+1))/2; %d(re(p))/dt=0 at final time
            p(2*((j+1)*NT-1)+2) = 0; %im(p)=0 at final time
        end
    else
        for j=0:(N-1) %changing phi to obey boundary conditions with new theta
            p(2*(j*NT+1)+1) = (p(2*j*NT+1) + p(2*(j*NT+1)+1))/2; %d(re(p))/dt=0 at initial time
            p(2*j*NT+1) = (p(2*j*NT+1) + p(2*(j*NT+1)+1))/2; %d(re(p))/dt=0 at initial time
            p(2*j*NT+2) = 0; %im(p)=0 at initial time
            p(2*((j+1)*NT-1)+1) = (p(2*((j+1)*NT-1)+1) + p(2*((j+1)*NT-2)+1))/2; %d(re(p))/dt=0 at final time
            p(2*((j+1)*NT-2)+1) = (p(2*((j+1)*NT-1)+1) + p(2*((j+1)*NT-2)+1))/2; %d(re(p))/dt=0 at final time
            p(2*((j+1)*NT-1)+2) = 0; %im(p)=0 at final time
        end
    end
        
    Xwait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (runsTest > closeness || runsCount<minRuns) && Xwait%beginning newton-raphson loop
        runsTest = abs(action - actionLast)/abs(actionLast); %note this won't work if action goes to zero
        runsCount = runsCount + 1;
        actionLast = action;
        
        minusDS = zeros(2*Tdim+1,1); %-dS/d(p)
        Cp = vecComplex(p,Tdim); %complex p, excluding real lagrange multiplier
        Chi0 = complex(zeros(N,1)); %to fix zero mode, alla kuznetsov, dl[7], though not quite
        %%CpNeg = complex(zeros(Tdim,1));

        for j=0:(N-1) %explicitly 2d, evaluating Chi0, equal to the zero mode at t=(Nt-1)
            if j~=(N-1)
                Chi0(j+1) = (p(2*((j+2)*Nt-1)+1)+1i*p(2*((j+2)*Nt-1)+1)-p(2*((j+1)*Nt-1)+1)+1i*p(2*((j+1)*Nt-1)+1))/a; %%note only use real derivative - this is a fudge due to initial input
            else
                Chi0(j+1) = (p(2*(Nt-1)+1)+1i*p(2*(Nt-1)+1)-p(2*(N*Nt-1)+1)+1i*p(2*(N*Nt-1)+1))/a;
            end
        end
        Chi0 = Chi0/norm(Chi0);
        
        nonz = 0; %number of nonzero elements of DDS
        if inP == 'b' || inP == 't' || inP == 'f'
            nonz = 5*N^(d-1) + 4*N^(d-1)*(Nt-2)*(2*d+1) + 4*Tdim;
        elseif inP == 'p'
            nonz = 6*N^(d-1) + 4*N^(d-1)*(Nt-2)*(2*d+1) + 4*Tdim;
        else %number of nonzero elements when fixing zero mode and with full boundary conditions - check
            nonz = 5*Tdim*N^(d-1) + N^(d-1) + 4*(Nt-2)*N^(d-1)*(2*d+1);
        end
        DDSm = zeros(nonz,1); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonz,1); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonz,1); %values of non-zero elements of DDS - don't forget to initialize DDS
        
        action = complex(0); %initializing to zero
        kinetic = complex(0);
        potL = complex(0);
        potE = complex(0);
        clear c3;
        
        for j=0:(Tdim-1)
            t = intCoord(j,0,Nt);
            x = intCoord(j,1,Nt);
            dtj = dt(j);
            siteMeasure = a*Dt(j); %for sites in time
            linkMeasure = a*dtj; %for links in time
            
            potL = potL - siteMeasure*(lambda/8)*(Cp(j+1)^2-v^2)^2;
            potE = potE - siteMeasure*epsilon*(Cp(j+1)-v)/v/2;
            if x~=(N-1) && x<(N-1)%evaluating spatial kinetic part
                kinetic = kinetic - siteMeasure*(Cp(j+Nt+1)-Cp(j+1))^2/a^2/2;
            elseif x~=(N-1) %avoinding using neigh and modulo as its slow
                kinetic = kinetic - siteMeasure*(Cp(j+1-Nt*(N-1))-Cp(j+1))^2/a^2/2;
            end
            if t==(Nt-1)
                if inP=='p' || inP=='b' %boundary conditions
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1/b; %zero time derivative %don't forget to clear c3
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j-1)+1; DDSv(c3) = -1/b;
                elseif inP=='t' || inP=='f'
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1; %zero change at boundary
                    DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1;
                end
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*Tdim+1; DDSv(c3) = real(siteMeasure*Chi0(x+1)); %zero mode lagrange constraint
                DDSm(c3) = 2*j+2; DDSn(c3) = 2*Tdim+1; DDSv(c3) = imag(siteMeasure*Chi0(x+1)); %the constraint is real but its derivative wrt phi may be complex
            else
                kinetic = kinetic + linkMeasure*((Cp(j+2) - Cp(j+1))/dtj)^2/2;
                if t==0
                    if inP=='b' || inP=='t' || inP=='f'
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1; %zero change at boundary
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1;
                    elseif inP=='p'
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = -1/b; %zero time derivative
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part zero
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+1)+1; DDSv(c3) = 1/b;
                    end
                else
                    dtjm = dt(j-1);
                    for k=0:(2*d-1)
                        sign = (-1)^k;
                        %%deltaSign = (sign-1)/2; %deltaSign=0 if sign=+1 and deltaSign=-1 if sign=-1
                        direc = floor(k/2);
                        if direc == 0
                            dtd = dtj;
                            if sign==-1
                                dtd = dtjm;
                            end
                            minusDS(2*j+1) = minusDS(2*j+1) + real(a*Cp(j+sign+1)/dtd);
                            minusDS(2*j+2) = minusDS(2*j+2) + imag(a*Cp(j+sign+1)/dtd);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+sign)+1; DDSv(c3) = -real(a/dtd);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+sign)+2; DDSv(c3) = imag(a/dtd);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*(j+sign)+1; DDSv(c3) = -imag(a/dtd);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*(j+sign)+2; DDSv(c3) = -real(a/dtd);
                        else
                            neighb = neigh(j,direc,sign,Nt);
                            minusDS(2*j+1) = minusDS(2*j+1) - real(siteMeasure*Cp(neighb+1)/a^2); %note Dt(j) = Dt(j+Nt) etc.
                            minusDS(2*j+2) = minusDS(2*j+2) - imag(siteMeasure*Cp(neighb+1)/a^2);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*neighb+1; DDSv(c3) = real(siteMeasure/a^2);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*neighb+2; DDSv(c3) = -imag(siteMeasure/a^2);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+1; DDSv(c3) = imag(siteMeasure/a^2);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+2; DDSv(c3) = real(siteMeasure/a^2);
                        end
                    end
                    temp0 = a*(1/dtj + 1/dtjm);
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
        for j=0:(N-1) %adding last two rows, with lagrange multiplier terms
            minusDS(2*Tdim+1) = minusDS(2*Tdim+1) - real(a*b*Chi0(j+1)*Cp((j+1)*Nt));
            
            DDSm(c3) = 2*Tdim+1; DDSn(c3) = 2*((j+1)*Nt-1)+1; DDSv(c3) = real(a*b*Chi0(j+1)); %at t=Nt-1
            DDSm(c3) = 2*Tdim+1; DDSn(c3) = 2*((j+1)*Nt-1)+2; DDSv(c3) = -imag(a*b*Chi0(j+1));                          
        end
        clear c3; %returning persistent output of c3 to 1
        kinetic = 2*kinetic; potL = 2*potL; potE = 2*potE; %as we only calculated half of bubble in time direction
        action = kinetic + potL + potE;
        DDSm(DDSv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
        DDSn(DDSv==0) = [];
        DDSv(DDSv==0) = [];
        DDS = sparse(DDSm,DDSn,DDSv,2*Tdim+1,2*Tdim+1);

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
                pause(5);
                [DDSm,DDSn,DDSv] = find(DDS);
                save data/DDS.mat DDSm DDSn DDSv;
                disp(['printed DDS in data/DDS.mat on run ',num2str(runsCount)]);
            else
                disp('early print error');
            end
        end

        [orderRow,orderCol,r,s,cc,rr] = dmperm(DDS); %preordering - gets vector order (and perhaps a second vector) - options are colamd, colperm and dmperm (which may produce 2 ordering vectors)
        
        setup.type = 'ilutp'; %preconditioning - incomplete LU factorization does not increase the number of non-zero elements in DDS - options are 'nofill', 'ilutp' and 'crout'
        setup.droptol = 1e-6; %drop tolerance is the minimum ratio of (off-diagonal) abs(U_ij) to norm(DDS(:j))
        setup.thresh = 0; %if 1 forces pivoting on diagonal, 0 to turn off
        setup.udiag = 0; %if 1 this replaces zeros in upper diagonal with droptol, 0 to turn off
        [Lo,Up] = ilu(DDS(orderRow,orderCol),setup);

        tol = 1e-6; %tolerance for norm(Ax-b)/norm(b), consider increasing if procedure is slow
        maxit = 50; %max number of iterations
        delta = zeros(2*Tdim+1);
        [delta(orderCol),flag,relres,iter,resvec] = lsqr(DDS(orderRow,orderCol),minusDS(orderRow),tol,maxit,Lo,Up); %finding solution iteratively. consider changing bicg to bicgstab, bicgstabl, cgs, gmres, lsqr, qmr or tfqmr 
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

        p = p + delta(orderCol)'; %p -> p'

        %pNeg and pZer plus log(det(DDS)) stuff

        stopTime = toc;
        [Xwait,aq] = convergenceQuestions(runsCount, runsTest, aq, stopTime, action); %discovering whether or not n-r has converged, and stopping if it is wildly out

    end %closing newton-raphson loop

    %propagating solution back in minkowskian time
    %checking that converged
    
    if loop==0 %printing to terminal
        fprintf('%12s','time', 'runs','d','N','X','re(action)','im(action)'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('\n');
    end
    fprintf('%12g',toc,runsCount,d,N,X,real(action));
    fprintf('%12g\n',imag(action));
    
    actionOut = fopen('data/picAction.dat','a'); %saving action etc to file
    fprintf(actionOut,'%12g',toc,runsCount,d,N,X,real(action));
    fprintf(actionOut,'%12g\n',imag(action));
    fclose(actionOut);
    
    save( ['/data/picVectors',num2str(loop),'.mat'], Cp, minusDS);%saving phi and minusDS to file
    
end%closing while loop