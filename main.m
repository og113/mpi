%script to solve boundary value problem with non-zero T and angle. input is
%the periodic instanton and then steps increasing angle or T.

global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim; %defining global variables
global R X lambda mass v epsilon angle theta;

minFileNo = input('which data/picOut#.mat file to load first? (#) '); %loading periodic instanton
maxFileNo = input('which data/picOut#.mat file to load last? (#) ');
for fileNo = minFileNo:maxFileNo
data = load(['data/picOut',num2str(fileNo),'.mat']);
data.DDS = []; data.minusDS = []; %freeing some memory

eigenData = load('data/eigens.mat'); %NEED TO CHANGE THIS ONCE WE CAN LOOP IN FILENO
negVal = eigenData.D;
negVec = eigenData.V;
%negVal = input('input vacuum bubble negative eigenvalue: ');
%negVec = input('input vacuum bubble negative eigenvector: ');

disp(['Lb = T/2 = ', num2str(data.Lb)]); %asking input questions
disp('angle = 0');
maxTheta = input('input final value of angle (input 0 for no steps) ');
%%maxLt = input('and input final value of Lb');
totalLoops = 1;
if maxTheta~=0
    totalLoops = input('and number of steps to get there ');
end

parametersMain(data); %assigning global variables according to data and parametersMain.m

Cp = data.tCp; data.tCp = []; %assigning phi from pic.m output, and freeing up memory
p = data.tp; data.tp = [];
p(end+1) = v; %second lagrange multiplier to remove time zero mode

Omega = omega(N); %for implementation of initial time boundary conditions
eOmega = Eomega(N); %comes with an extra power of the energy

for loop=0:(totalLoops-1) %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if totalLoops>1
        theta = loop*maxTheta/(totalLoops - 1);
    end
    gamma = exp(-theta);
    
    tic; %starting the clock
    
    action = complex(2);
    
    actionLast = complex(1); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    actionTest = 1;
    vectorTest = 1;
    closenessA = 1e-4;
    closenessV = 1e-3;
    minRuns = 3;
        
    Xwait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (actionTest(end) > closenessA || vectorTest(end) > closenessV || runsCount<minRuns) && Xwait%beginning newton-raphson loop
        runsCount = runsCount + 1;
        
        minusDS = zeros(2*Tdim+2,1); %-dS/d(p)
        chiX = zeros(Nb*N,1); %to fix x zero mode, field made perpendicular to it
        chiT = zeros(2*NT*N+2,1); %to fix (real) t zero mode, time derivative of field made perpendicular to it

        for j=0:(N-1) %evaluating chiX, equal to the zero mode at t=(Nb-1)
            posB = (j+1)*Nb-1; %position on section b of contour
            posBf = (Na+Nb-1)+j*NT; %position of b on full contour
            posT = (j+1)*NT-1;
            chiX(posB+1) = p(2*neigh(posBf,1,1,NT)+1)-p(2*neigh(posBf,1,-1,NT)+1);
            %chiX(posB-1+1) = p(2*neigh(posBf-1,1,1,NT)+1)-p(2*neigh(posBf-1,1,-1,NT)+1);
            chiT(2*posT+1) = negVec(2*posB+1); %arbitrary atempt to fix zero mode
        end
        %chiX = v*chiX/norm(chiX);
        
        nonz = (2*2*2+2)*NT*(N-2) + 6*Tdim + 2*2*N +4*4*N; %number of nonzero elements of DDS - an overestimate
        DDSm = zeros(nonz,1); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonz,1); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonz,1); %values of non-zero elements of DDS - don't forget to initialize DDS
        
        kinetic = complex(0);%initializing to zero
        potL = complex(0);
        potE = complex(0);
        clear c3;
        
        for j=0:(Tdim-1)
            t = intCoord(j,0,NT);
            x = intCoord(j,1,NT);
            dtj = tdt(j);
            siteMeasure = a*tDt(j); %for sites in time
            linkMeasure = a*dtj; %for links in time
            
            if (t>=Na) && (t<(Na+Nb)) %fixing zero modes
                posB = t-Na+x*Nb;
                if abs(chiX(posB+1))>1e-16 %orthogonal to p
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*Tdim+1; DDSv(c3) = a*chiX(posB+1);
                    DDSm(c3) = 2*Tdim+1; DDSn(c3) = 2*j+1; DDSv(c3) = a*chiX(posB+1);  
                    minusDS(2*j+1) = minusDS(2*j+1) - a*chiX(posB+1)*p(2*Tdim+1);
                    minusDS(2*Tdim+1) = minusDS(2*Tdim+1) - a*chiX(posB+1)*p(2*j+1);
                end
            end
            if abs(chiT(2*j+1))>1e-16 %orthogonal to dp/dt
                DDSm(c3) = 2*(j+1)+1; DDSn(c3) = 2*Tdim+2; DDSv(c3) = a*chiT(2*j+1); %there should be no chiT on the final time slice or this line will go wrong
                DDSm(c3) = 2*Tdim+2; DDSn(c3) = 2*(j+1)+1; DDSv(c3) = a*chiT(2*j+1);
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*Tdim+2; DDSv(c3) = -a*chiT(2*j+1);
                DDSm(c3) = 2*Tdim+2; DDSn(c3) = 2*j+1; DDSv(c3) = -a*chiT(2*j+1);
                minusDS(2*(j+1)+1) = minusDS(2*(j+1)+1) - a*chiT(2*j+1)*p(2*Tdim+2);
                minusDS(2*Tdim+2) = minusDS(2*Tdim+2) - a*chiT(2*j+1)*p(2*j+1);
                minusDS(2*j+1) = minusDS(2*j+1) + a*chiT(2*j+1)*p(2*Tdim+2);
                minusDS(2*Tdim+2) = minusDS(2*Tdim+2) + a*chiT(2*j+1)*p(2*j+1);
            end
            
            potL = potL - siteMeasure*(lambda/8)*(Cp(j+1)^2-v^2)^2;
            potE = potE - siteMeasure*epsilon*(Cp(j+1)-v)/v/2;
            kinetic = kinetic - siteMeasure*(Cp(neigh(j,1,1,NT)+1)-Cp(j+1))^2/a^2/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% boundary conditions
            if t==(NT-1)
                %zero eigenvalue due to these boundary conditions
                DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part of field at final time boundary
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j-1)+2; DDSv(c3) = 1; %zero imaginary part of time derivative at final time boundary - with other condition               z
            else
                kinetic= kinetic + linkMeasure*((Cp(j+2) - Cp(j+1))/dtj)^2/2;
                if t==0
                    if abs(theta)<1e-16
                        %zero eigenvalue due to these boundary conditions
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+1)+2; DDSv(c3) = 1; %zero imaginary time derivative - with below condition
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part
                    else
                        for k=1:(2*d-1)
                            sign = (-1)^(k-1);
                            %%deltaSign = (sign-1)/2; %deltaSign=0 if sign=+1 and deltaSign=-1 if sign=-1
                            direc = floor(k/2);
                            if direc == 0
                                minusDS(2*j+1) = minusDS(2*j+1) + real(a*Cp(j+sign+1)/dtj);
                                minusDS(2*j+2) = minusDS(2*j+2) + imag(a*Cp(j+sign+1)/dtj);
                            else
                                neighb = neigh(j,direc,sign,NT);
                                minusDS(2*j+1) = minusDS(2*j+1) - real(siteMeasure*Cp(neighb+1)/a^2); %note Dt(j) = Dt(j+Nt) etc.
                                minusDS(2*j+2) = minusDS(2*j+2) - imag(siteMeasure*Cp(neighb+1)/a^2);
                            end
                        end
                        temp0 = a/dtj;
                        temp1 = siteMeasure*(2*(d-1)*Cp(j+1)/a^2 + (lambda/2)*Cp(j+1)*(Cp(j+1)^2-v^2) + epsilon/2/v);

                        minusDS(2*j+1) = minusDS(2*j+1) + real(temp1 - temp0*Cp(j+1));
                        minusDS(2*j+2) = minusDS(2*j+2) + imag(temp1 - temp0*Cp(j+1));

                        for k=0:(N-1)
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*k*NT+1; DDSv(c3) = real(1i*Omega(x+1,k+1)*(1+gamma)/(1-gamma));
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*k*NT+2; DDSv(c3) = real(-Omega(x+1,k+1)*(1-gamma)/(1+gamma));
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*k*NT+1; DDSv(c3) = imag(1i*Omega(x+1,k+1)*(1+gamma)/(1-gamma));
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*k*NT+2; DDSv(c3) = imag(-Omega(x+1,k+1)*(1-gamma)/(1+gamma));
                        end
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of boundary conditions
                else
                    dtjm = tdt(j-1);
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
                            neighb = neigh(j,direc,sign,NT);
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
        clear c3; %returning persistent output of c3 to 1
        action = kinetic + potL + potE;
        DDSm(DDSv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
        DDSn(DDSv==0) = [];
        DDSv(DDSv==0) = [];
        DDS = sparse(DDSm,DDSn,DDSv,2*Tdim+2,2*Tdim+2);
        
        undetermined = length(DDS)-sprank(DDS);
        if undetermined>1e-16
           [chiT,D] = eigs(DDS,1,1e-10);
           continue %goes to beginning of while loop and starts again
        end
        
        erg = 0;
        num = 0;
        W = 0;
        if abs(theta)<1e-16
            for j=0:(N-1)
                for k=0:(N-1)
                    num = num + Omega(j+1,k+1)*p(2*j*NT+1)*p(2*k*NT+1) ...
                        - Omega(j+1,k+1)*p(2*j*NT+2)*p(2*k*NT+2); %the sign here may be wrong - i feel intuitively that it aught to be positive
                    erg = erg + eOmega(j+1,k+1)*p(2*j*NT+1)*p(2*k*NT+1) ...
                        - eOmega(j+1,k+1)*p(2*j*NT+2)*p(2*k*NT+2); %likewise with sign
%extra boundary term vanishes for the periodic instanton (theta = 0)
                end
            end
        else
            for j=0:(N-1)
                for k=0:(N-1)
                    num = num + 2*gamma*Omega(j+1,k+1)*p(2*j*NT+1)*p(2*k*NT+1)/(1+gamma)^2 ...
                        + 2*gamma*Omega(j+1,k+1)*p(2*j*NT+2)*p(2*k*NT+2)/(1-gamma)^2;
                    erg = erg + 2*gamma*eOmega(j+1,k+1)*p(2*j*NT+1)*p(2*k*NT+1)/(1+gamma)^2 ...
                        + 2*gamma*eOmega(j+1,k+1)*p(2*j*NT+2)*p(2*k*NT+2)/(1-gamma)^2;
                    W = W - (1-gamma)*Omega(j+1,k+1)*p(2*j*NT+1)*p(2*k*NT+1)/(1+gamma) ...
                        + (1+gamma)*Omega(j+1,k+1)*p(2*j*NT+2)*p(2*k*NT+2)/(1-gamma);
                end
            end
        end
        W = W + num*theta + erg*2*Lb - 2*imag(action);
        W = -lambda*W;
        
        
        small = norm(minusDS); %normalising problem
        smaller = small/(2*Tdim+2);
        cutoff = 1e-16;
        normed = 0;
        if smaller < cutoff
           Xwait = 0;
           break;
        else
            normed = 1;
            minusDS = minusDS/small;
            DDS = DDS/small;
        end
        
        limit = 1e12;
        c = condest(DDS); %finding estimate and bound for condition number
        eigen = eigs(DDS,1,1e-6);
        if c>limit || abs(eigen)<1/limit || undetermined>1e-16;
            singular = svds(DDS,1,0);
            fprintf('%12s','size - structural rank = ');
            fprintf('%12g\n',undetermined);
            fprintf('%12s','condest = ');
            fprintf('%12g\n',c);
            fprintf('%12s','smallest singular value is = ');
            fprintf('%12g\n',singular);
            fprintf('%12s','smallest eigenvalue is = ');           
            fprintf('%12g\n',eigen);
        end
        
        %[orderRow,orderCol,r,s,cc,rr] = dmperm(DDS); %preordering - gets vector order (and perhaps a second vector) - options are colamd, colperm and dmperm (which may produce 2 ordering vectors)
        
        %setup.type = 'ilutp'; %preconditioning - incomplete LU factorization does not increase the number of non-zero elements in DDS - options are 'nofill', 'ilutp' and 'crout'
        %setup.droptol = 1e-6; %drop tolerance is the minimum ratio of (off-diagonal) abs(U_ij) to norm(DDS(:j))
        %setup.thresh = 0; %if 1 forces pivoting on diagonal, 0 to turn off
        %setup.udiag = 0; %if 1 this replaces zeros in upper diagonal with droptol, 0 to turn off
        %[Lo,Up] = ilu(DDS(orderRow,orderCol),setup);

        %tol = 1e-6; %tolerance for norm(Ax-b)/norm(b), consider increasing if procedure is slow
        %maxit = 50; %max number of iterations
        %[delta(orderCol),flag,relres,iter,resvec] = lsqr(DDS(orderRow,orderCol),minusDS(orderRow),tol,maxit,Lo,Up); %finding solution iteratively. consider changing bicg to bicgstab, bicgstabl, cgs, gmres, lsqr, qmr or tfqmr 
        flag = 0;
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

        delta = DDS\minusDS; %this is where the magic happens - direct gaussian elimination
        
        if normed
            minusDS = minusDS*small;
            DDS = DDS*small;
        end
        
        vectorTest = [vectorTest, norm(delta)/norm(p)];
        actionTest = [actionTest, abs(action - actionLast)/abs(actionLast)]; %note this won't work if action goes to zero
        actionLast = action;
          
        if size(delta,1)==1
            p = p + delta'; %p -> p'
        else
            p = p + delta;
        end
        
        Cp = vecComplex(p,Tdim); 

        stopTime = toc;
        [Xwait] = convergenceQuestions(runsCount, stopTime, action); %discovering whether or not n-r has converged, and stopping if it is wildly out

        save( ['data/mainEarly',num2str(loop),num2str(runsCount),'.mat'], 'p', 'DDS', 'Cp', 'minusDS', 'd', 'N', 'Na', 'Nb', 'Nc', 'NT', 'lambda', 'mass', 'R', 'Lb','L');
        if 1 == 1 %loop==0
            disp(['runscount : ',num2str(runsCount),', time: ',num2str(toc),', actionTest: ',num2str(actionTest(end)),', vectorTest: ',num2str(vectorTest(end))]); %just to see where we are for first run
        end
        
    end %closing newton-raphson loop
    
    if 1==1 %loop==0 %printing to terminal
        fprintf('%6s','time', 'runs','N','Na','Nb', 'Nc', 'mass','lambda','R','Lb','theta'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('%12s','num','erg','re(action)','im(action)','W');
        fprintf('\n');
    end
    fprintf('%6g',toc,runsCount,N,Na,Nb,Nc,mass,lambda,R,Lb,theta);
    fprintf('%12g',num,erg,real(action),imag(action),W);
    fprintf('\n');
    
    actionOut = fopen('data/picAction.dat','a'); %saving action etc to file
    fprintf(actionOut,'%12g',toc,runsCount,N,X,Lb,theta,num,erg,real(action),imag(action),W);
    fprintf(actionOut,'\n');
    fclose(actionOut);
    
    save( ['data/main',num2str(fileNo),'_',num2str(loop),'.mat'], 'p', 'DDS', 'Cp', 'minusDS', 'd', 'N', 'Na', 'Nb', 'Nc', 'NT', 'lambda', 'mass', 'R','Lb','L','theta');%saving phi and minusDS to file
    
end%closing parameter loop
end %closing fileNo loop

data = load(['data/main',num2str(fileNo),'_',num2str(loop),'.mat']);
data.tCp = data.Cp;
plotTphi(data);