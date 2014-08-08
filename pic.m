%script to find the periodic instanton of the stable phi^4 with negative
%mass

global d N Na Nb Nc L Ltemp La Lb Lc a b Adim Bdim Cdim Tdim; %defining global variables
global R X lambda mass v epsilon angle amp;

aq.inputP = 'i'; %struct to hold answers to questions aq short for 'answers to questions' - defauLbs in initialization
aq.perturbResponse = 'n';
aq.loopResponse = 'n';
aq.parameterChoice = 'N';
aq.minValue = 32;
aq.maxValue = 64;
aq.totalLoops = 5;

inP = aq.inputP; %just because I write this a lot
parameters(inP); %assigning global variables according to parameters.m

fprintf('%10s','N', 'Na','Nb','Nc','L','Lb','R','mass','lambda','epsilon'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
fprintf('\n');
fprintf('%10g',N,Na,Nb,Nc,L,Lb,R,mass,lambda);
fprintf('%10g\n',epsilon);

tempAq = askQuestions; %asking questions
fields = fieldnames(tempAq);
for j=1:numel(fields) %using non-empty responses to questions to fill in aq
    if ~isempty(tempAq.(fields{j}))
       aq.(fields{j}) = tempAq.(fields{j}); 
    end
end
inP = aq.inputP; %just because I write this a lot
parameters(inP);

negVal = 0;
negVec = zeros(2*N*Nb+1,1);
if strcmp(inP,'q') || strcmp(inP,'p') || strcmp(inP,'i')
    eigenData = load('data/eigens.mat');
    negVal = eigenData.D;
    negVec = eigenData.V;
    %negVal = input('input vacuum bubble negative eigenvalue: ');
    %negVec = input('input vacuum bubble negative eigenvector: ');
end
if strcmp(inP,'i')
    inputData = load('data/picOut0.mat');
end

for loop=0:(aq.totalLoops-1) %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if strcmp(aq.loopResponse,'y') && strcmp(aq.parameterChoice,'N')
        N = aq.minValue + floor(loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1));
        changeParameters(N,'N',inP);
    elseif strcmp(aq.loopResponse,'y')
        loopParameter = aq.minValue + loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1);
        changeParameters (loopParameter,aq.parameterChoice, inP);
    end
    if strcmp(inP,'i') && loop>0
        inputData = load(['data/picOut',num2str(loop-1),'.mat']);
    end
    tic; %starting the clock
    
    S1 = 2*mass^3/3/lambda; %this is twice the value in the coleman paper
    twAction = -solidAngle(d)*epsilon*R^d/d + solidAngle(d)*R^(d-1)*S1; %thin-wall bubble action
    alpha = 10; %determines range over which tanh(x) is used
    if ~strcmp(aq.parameterChoice,'amp')
        amp = 2*(Lb-R)/R; %determines admixture of negative mode - trial and error
    end
    action = complex(2);
    
    actionLast = complex(twAction/2); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    actionTest = 1;
    deltaTest = 1;
    normTest = 1;
    maxTest = 1;
    closenessA = 1;
    closenessD = 1;
    closenessN = 1e-5;
    closenessM = 1e-4;
    minRuns = 3;
    
    p = zeros(2*Bdim+1,1); %phi, in the euclidean domain
    perturbReal = zeros(Bdim,1); %perturbations
    perturbImag = zeros(Bdim,1);
    %%pNeg = zeros(2*Bdim,1); %negative eigenvector
    
    %syms x
    %minima = vpasolve(x^3 -v^2*x + epsilon/v/lambda,x); %solving V'(p)=0
    polynomial = [1, 0 , -v^2, epsilon/v/lambda];
    minima = roots(polynomial);
    minima = sort(minima); %roots sorted in ascending order https://www.google.co.uk/search?client=ubuntu&channel=fs&q=matlab+random+square+that+i+can%27t+click+in&ie=utf-8&oe=utf-8&gl=uk&gws_rd=cr&ei=L8dxU8LPCo_d7Qbg3IGYCg#channel=fs&gl=uk&q=matlab+annoying+square+that+i+can%27t+click+in
    
    if ~strcmp(aq.perturbResponse,'n') %assigning values to perturbations if user wants perturbations
        perturbReal = v*1e-4*rand(Bdim,1);
        perturbImag = v*1e-4*rand(Bdim,1);
        for j=1:(d-1)
            for k=0:(N-1)
                if inP == 't' || inP =='f' %pure vacua
                    perturbReal(k*Nb*N^(j-1)+1) = 0; %zero perturbation at both (euclidean) time boundaries
                    perturbReal((k+1)*Nb*N^(j-1)) = 0;
                    perturbImag(k*Nb*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nb*N^(j-1)) = 0;
                elseif inP == 'b' %bubble
                    perturbReal(k*Nb*N^(j-1)+1) = 0; %zero perturbation at initial time
                    perturbReal((k+1)*Nb*N^(j-1)) = perturbReal((k+1)*Nb*N^(j-1)-1); %zero derivative at final time
                    perturbImag(k*Nb*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nb*N^(j-1)) = 0;
                elseif inP == 'p' || inP == 'q' || inP == 'i' %periodic instanton - 'q' is the other approximation to the periodic instanton
                    perturbReal(k*Nb*N^(j-1)+1) = perturbReal(k*Nb*N^(j-1)+2); %zero time derivative at both time boundaries
                    perturbReal((k+1)*Nb*N^(j-1)) = perturbReal((k+1)*Nb*N^(j-1)-1);
                    perturbImag(k*Nb*N^(j-1)+1) = 0;
                    perturbImag((k+1)*Nb*N^(j-1)) = 0;
                end
            end
        end
    end
    
    for j=0:(Bdim-1) %assigning input phi according to inputP, note that the periodic instanton input is explicitly 2d
        t = eCoord(j,0);
        x = eCoord(j,1);
        rhoSqrd = -t^2;
        rho1Sqrd = -t^2;
        rho2Sqrd = -t^2; 
        rhoSqrd = rhoSqrd + x^2;
        rho1Sqrd = rho1Sqrd + (x+R*cos(angle))^2;
        rho2Sqrd = rho2Sqrd + (x-R*cos(angle))^2;
        rho = real(sqrt(rhoSqrd)); %rho should be real even without the real()
        rho1 = real(sqrt(rho1Sqrd));
        rho2 = real(sqrt(rho2Sqrd));
        if R<alpha/mass
            disp(['X = R*mass is too small. not possible to give thinwall input. it should be less that ',num2str(alpha)]);
        else
            if strcmp(inP,'t')
                p(2*j+1) = minima(1);
            elseif strcmp(inP,'f')
                p(2*j+1) = minima(3);
            elseif strcmp(inP,'b') || strcmp(inP,'q')
                if rho<(R-alpha/mass)
                    p(2*j+1) = minima(1);
                elseif rho>(R+alpha/mass)
                    p(2*j+1) = minima(3);
                else
                    p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(3)-minima(1))*tanh(mass*(rho-R)/2)/2.0;
                    %%pNeg(2*j+1) = v/cosh(mass*(rho-R)/2)^2;
                end
            elseif strcmp(inP,'p')
                if rho1<(R-alpha/mass) && rho2<(R-alpha/mass)
                    p(2*j+1) = minima(1);
                elseif rho1>(R+alpha/mass) || rho2>(R+alpha/mass)
                    p(2*j+1) = minima(3);
                elseif real(eCoord(j,1))>0 %explicitly 2d here, note that the coord should be real
                    p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(3)-minima(1))*tanh(mass*(rho1-R)/2)/2.0;
                    %%pNeg(2*j+1) = v/cosh(mass*(rho1-R)/2)^2;
                elseif real(eCoord(j,1))<0
                    p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(3)-minima(1))*tanh(mass*(rho2-R)/2)/2.0;
                    %%pNeg(2*j+1) = v/cosh(mass*(rho2-R)/2)^2;
                else
                    p(2*j+1) = minima(1); %if eCoord(j,1) == 0
                end
            elseif strcmp(inP,'i')
                p(2*j+1) = inputData.p(2*j+1);
            end
            if strcmp(inP,'q') || (strcmp(inP,'i') && strcmp(inputData.inP,'b'))
                p(2*j+1) = p(2*j+1) + amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+1); %adding the negative eigenvector
                p(2*j+2) = p(2*j+2) + amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+2);
            end
            if aq.perturbResponse == 'r' || aq.perturbResponse =='b'
                p(2*j+1) = p(2*j+1) + perturbReal(j+1);
            end
            if aq.perturbResponse == 'i' || aq.perturbResponse =='b'
                p(2*j+2) = p(2*j+2) + perturbImag(j+1);
            end
        end
    end
    
    p(Bdim+1) = v; %initializing Lagrange parameter for dp/dx zero mode
    
    if (strcmp(inP,'p') || strcmp(inP,'q')) && 1 %fixing input periodic instanton to have zero time derivative at time boundaries
        open = 0; %value of 0 assigns all weight to boundary, value of 1 to neighbour of boundary
        for j=0:(N-1)
            p(2*j*Nb+1) = (1.0-open)*p(2*j*Nb+1) + open*p(2*(j*Nb+1)+1);%intiial time real
            p(2*(j*Nb+1)+1) = p(2*j*Nb+1);
            p(2*j*Nb+2) = (1.0-open)*p(2*j*Nb+2) + open*p(2*(j*Nb+1)+2);%initial time imag
            p(2*(j*Nb+1)+2) = p(2*j*Nb+2);
            p(2*((j+1)*Nb-1)+1) = open*p(2*((j+1)*Nb-2)+1) + (1.0-open)*p(2*((j+1)*Nb-1)+1);%final time real
            p(2*((j+1)*Nb-2)+1) = p(2*((j+1)*Nb-1)+1);
            p(2*((j+1)*Nb-2)+2) = open*p(2*((j+1)*Nb-1)+2) + (1.0-open)*p(2*((j+1)*Nb-2)+2);%final time imag
            p(2*((j+1)*Nb-1)+2) = p(2*((j+1)*Nb-2)+2);
        end
    end
    
    %%if inP == 'b' || inP == 'p' %fixing norm of pNeg
        %%pNeg = pNeg/norm(pNeg);
    %%end
        
    Xwait = 1; %this is just a parameter to stop when things are slow perhaps because the newton-raphson loop isn't converging
    while (actionTest(end) > closenessA || deltaTest(end) > closenessD || normTest(end) > closenessN || maxTest(end) > closenessM || runsCount<minRuns) && Xwait%beginning newton-raphson loop
        runsCount = runsCount + 1;
        
        minusDS = zeros(2*Bdim+1,1); %-dS/d(p)
        Cp = vecComplex(p,Bdim); %complex p, excluding real lagrange muLbiplier
        Chi0 = zeros(N*Nb,1); %to fix zero mode, alla kuznetsov, dl[7], though not quite
        %%CpNeg = complex(zeros(Bdim,1));

        for j=0:(N-1) %explicitly 2d, evaluating Chi0, equal to the zero mode at t=(Nb-1)
            pos = (j+1)*Nb-1;
            %Chi0(pos+1) = negVec(2*pos+1);
            Chi0(pos+1) = p(2*neigh(pos,1,1,Nb)+1)-p(2*neigh(pos,1,-1,Nb)+1);
            %Chi0(pos-1+1) = p(2*neigh(pos-1,1,1,Nb)+1)-p(2*neigh(pos-1,1,-1,Nb)+1);
        end
        %Chi0 = v*Chi0/norm(Chi0);
        
        nonz = 0; %number of nonzero elements of DDS
        if inP == 'b' || inP == 't' || inP == 'f'
            nonz = 5*N^(d-1) + 4*N^(d-1)*(Nb-2)*(2*d+1) + 4*Bdim;
        elseif inP == 'p' || inP == 'q' || inP == 'i'
            nonz = 6*N^(d-1) + 4*N^(d-1)*(Nb-2)*(2*d+1) + 4*Bdim;
        else %number of nonzero elements when fixing zero mode and with full boundary conditions - check
            nonz = 5*Bdim*N^(d-1) + N^(d-1) + 4*(Nb-2)*N^(d-1)*(2*d+1);
        end
        DDSm = zeros(nonz,1); %row numbers of non-zero elements of DDS
        DDSn = zeros(nonz,1); %column numbers of non-zero elements of DDS
        DDSv = zeros(nonz,1); %values of non-zero elements of DDS - don't forget to initialize DDS
         
        kinetic = complex(0); %initializing to zero
        potL = complex(0);
        potE = complex(0);
        clear c3;
        
        for j=0:(Bdim-1)
            t = intCoord(j,0,Nb);
            x = intCoord(j,1,Nb);
            dtj = dt(j);
            Dtj = Dt1(j);
            siteMeasure = a*Dtj; %for sites in time
            linkMeasure = a*dtj; %for links in time
            
            if abs(Chi0(j+1))>1e-16 %zero mode lagrange constraint
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*Bdim+1; DDSv(c3) = a*Chi0(j+1);
                DDSm(c3) = 2*Bdim+1; DDSn(c3) = 2*j+1; DDSv(c3) = a*Chi0(j+1);  
                minusDS(2*j+1) = minusDS(2*j+1) - a*Chi0(j+1)*p(2*Bdim+1);
                minusDS(2*Bdim+1) = minusDS(2*Bdim+1) - a*Chi0(j+1)*p(2*j+1);                      
            end
            
            potL = potL - siteMeasure*(lambda/8)*(Cp(j+1)^2-v^2)^2;
            potE = potE - siteMeasure*epsilon*(Cp(j+1)-v)/v/2;
            if x~=(N-1) && x<(N-1)%evaluating spatial kinetic part
                kinetic = kinetic - siteMeasure*(Cp(j+Nb+1)-Cp(j+1))^2/a^2/2;
            elseif x==(N-1) %avoinding using neigh and modulo as its slow
                kinetic = kinetic - siteMeasure*(Cp(j+1-Nb*(N-1))-Cp(j+1))^2/a^2/2;
            end
            if t==(Nb-1)
                if inP=='p' || inP=='q' || inP=='b' || inP == 'i'%boundary conditions
                    DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = 1/b; %zero time derivative
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
                    elseif inP=='p' || inP=='q' || inP =='i'
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*j+1; DDSv(c3) = -1/b; %zero time derivative
                        DDSm(c3) = 2*j+2; DDSn(c3) = 2*j+2; DDSv(c3) = 1; %zero imaginary part
                        DDSm(c3) = 2*j+1; DDSn(c3) = 2*(j+1)+1; DDSv(c3) = 1/b;
                    end
                else
                    dtjm = dt(j-1);
                    for k=0:(2*d-1)
                        sign = (-1)^k;
                        dtd = dtj;
                        if sign==-1
                           dtd = dtjm; 
                        end
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
                            neighb = neigh(j,direc,sign,Nb);
                            minusDS(2*j+1) = minusDS(2*j+1) - real(Dtj*Cp(neighb+1)/a); %note Dt1(j) = Dt1(j+Nb) etc.
                            minusDS(2*j+2) = minusDS(2*j+2) - imag(Dtj*Cp(neighb+1)/a);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*neighb+1; DDSv(c3) = real(Dtj/a);
                            DDSm(c3) = 2*j+1; DDSn(c3) = 2*neighb+2; DDSv(c3) = -imag(Dtj/a);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+1; DDSv(c3) = imag(Dtj/a);
                            DDSm(c3) = 2*j+2; DDSn(c3) = 2*neighb+2; DDSv(c3) = real(Dtj/a);
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
        DDS = sparse(DDSm,DDSn,DDSv,2*Bdim+1,2*Bdim+1);
        
        small = norm(minusDS); %normalising problem
        smaller = small/(2*Bdim+1);
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
        
        undetermined = length(DDS)-sprank(DDS);
        if undetermined>1e-16
            fprintf('%12s','size - structural rank = ');
            break %goes to beginning of while loop and starts again
        end
        
        limit = 1e15*min([closenessA,closenessD,closenessN,closenessM]); %if c>limit numerical errors will be significant
        c = condest(DDS); %finding estimate and bound for condition number
        if c>limit
            fprintf('%12s','condest = ');
            fprintf('%12g\n',c);
            %singular = svds(DDS,1,1e-10);
            %fprintf('%12s','smallest singular value is = ');
            %fprintf('%12g\n',singular);
        end
        

        %[orderRow,orderCol,r,s,cc,rr] = dmperm(DDS); %preordering - gets vector order (and perhaps a second vector) - options are colamd, colperm and dmperm (which may produce 2 ordering vectors)
        %orderRow = dmperm(DDS);
        %orderRow(end) = length(orderRow); %as it sometimes gives 0 at the end
        
        %setup.type = 'nofill';%'ilutp'; %preconditioning - incomplete LU factorization does not increase the number of non-zero elements in DDS - options are 'nofill', 'ilutp' and 'crout'
        %setup.droptol = 1e-6; %drop tolerance is the minimum ratio of (off-diagonal) abs(U_ij) to norm(DDS(:j))
        %setup.thresh = 0; %if 1 forces pivoting on diagonal, 0 to turn off
        %setup.udiag = 0; %if 1 this replaces zeros in upper diagonal with droptol, 0 to turn off
        %[Lo,Up] = ilu(DDS(orderRow,orderCol),setup);
        %[Lo,Up] = ilu(DDS(orderRow,:),setup);
        
        %tol = 1e-6; %tolerance for norm(Ax-b)/norm(b), consider increasing if procedure is slow
        %maxit = 50; %max number of iterations
        %x0 = rand(2*Bdim+1,1);
        %delta = zeros(2*Bdim+1,1);
        %[delta(orderCol),flag,relres,iter,resvec] = lsqr(DDS(orderRow,orderCol),minusDS(orderRow),tol,maxit,Lo,Up,x0); %finding solution iteratively. consider changing bicg to bicgstab, bicgstabl, cgs, gmres, lsqr, qmr or tfqmr 
        %[delta,flag,relres,iter,resvec] = lsqr(DDS(orderRow,:),minusDS(orderRow),tol,maxit,Lo,Up,x0);
        flag=0;
        if flag ~=0 %flag = 0 means bicg has converged, if getting non-zero flags, output and plot relres ([deLba,flag] -> [deLba,flag,relres])
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
        %fprintf('%12s','iter = ');
        %fprintf('%12g\n',iter);
        %fprintf('%12s','relres = ');
        %fprintf('%12g\n',relres);
        
        delta = DDS\minusDS; %this is where the magic happens - direct gaussian elimination
        
        if normed
            minusDS = minusDS*small;
            DDS = DDS*small;
        end
        
        deltaTest = [deltaTest, norm(delta)/norm(p)];
        actionTest = [actionTest, abs(action - actionLast)/abs(actionLast)]; %note this won't work if action goes to zero
        actionLast = action;
        normTest = [normTest, norm(minusDS)/norm(p)];
        maxTest = [maxTest, max(abs(minusDS))/max(abs(p))];
        
        if size(delta,1)==1
            p = p + delta'; %p -> p'
        else
            p = p + delta;
        end

        %pNeg and pZer plus log(det(DDS)) stuff

        stopTime = toc;
        [Xwait, Xaq] = convergenceQuestions(aq, runsCount, stopTime, action); %discovering whether or not n-r has converged, and stopping if it is wildly out
        aq = Xaq;
        
        save( ['data/picEarly',num2str(loop),num2str(runsCount),'.mat'],'p', 'Cp', 'minusDS','DDS','action', 'd', 'N', 'Na', 'Nb' , 'Nc', 'lambda', 'mass', 'R','L','La','Lb','Lc','inP');
        %action
        if 1 == 1 %loop==0
            disp(['runscount : ',num2str(runsCount),', time: ',num2str(toc),', actionTest: ',num2str(actionTest(end)),', deltaTest: ',num2str(deltaTest(end))...
                ,', normTest: ',num2str(normTest(end)),', maxTest: ',num2str(maxTest(end))]); %just to see where we are for first run
        end
        
    end %closing newton-raphson loop

    %propagating solution back in minkowskian time
    %A1. initialize mp==mphi using last point of ephi and zeros- use complex phi
    ap = complex(zeros(Adim+N,1)); %extra N is to include boundary point
    for j=0:(N-1)
        ap(j*(Na+1)+1) = p(2*(j*Nb)+1) + 1i*p(2*(j*Nb)+2);%%roots(1);%%
    end

    %A2. initialize vel - defined at half steps, first step being at t=-1/2,
    %vel(t+1/2) := (p(t+1)-p(t))/dt
    vel = complex(zeros(Adim+N,1));
    dtau = -b;
    Dt0 = dtau;%%b/2*(-1+1i*up); - this is surely wrong!!
    for j=0:(N-1)
        k = j*(Na+1);
        vel(k+1) = 0; %%due to boundary condition
    end

    
    %A3. initialize acc using phi and expression from equation of motion and zeros-
    %complex
    acc = complex(zeros(Adim+N,1));
    for j=0:(N-1)
        k = j*(Na+1);
        acc(k+1) = ((Dt0/a^2)*(ap(neigh(k,1,1,Na+1)+1)+ap(neigh(k,1,-1,Na+1)+1)-2*ap(k+1)) ...
            -(lambda*Dt0/2)*ap(k+1)*(ap(k+1)^2-v^2) - epsilon*Dt0/2/v)/dtau;
    end

    %A7. run loop
    for j=1:Na
        for k=0:(N-1)
            l = j+k*(Na+1);
            vel(l+1) = vel(l) + dtau*acc(l);
            ap(l+1) = ap(l) + dtau*vel(l+1);%
        end
        for k=0:(N-1)
            l = j+k*(Na+1);
            acc(l+1) = (1/a^2)*(ap(neigh(l,1,1,Na+1)+1)+ap(neigh(l,1,-1,Na+1)+1)-2*ap(l+1)) ...
            -(lambda/2)*ap(l+1)*(ap(l+1)^2-v^2) - epsilon/2/v;    
        end
    end

    %now along c
    %C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
    cp = complex(zeros(Cdim+N,1)); %extra N is to include corner point
    for j=0:(N-1)
        cp(j*(Nc+1)+1) = p(2*(j*Nb+Nb-1)+1) + 1i*p(2*(j*Nb+Nb-1)+2);%%roots(1);%%
    end

    %C3. initialize vel - defined at half steps, first step being at t=-1/2,
    %vel(t+1/2) := (p(t+1)-p(t))/dt
    vel = complex(zeros(Cdim+N,1));
    dtau = b;
    Dt0 = dtau;%%b/2*(-1+1i*up); - this is surely wrong!!
    for j=0:(N-1)
        k = j*(Nc+1);
        vel(k+1) = 0; %%due to boundary condition
    end

    %C4. initialize acc using phi and expression from equation of motion and zeros-
    %complex
    acc = complex(zeros(Cdim+N,1));
    for j=0:(N-1)
        k = j*(Nc+1);
        acc(k+1) = ((Dt0/a^2)*(cp(neigh(k,1,1,Nc+1)+1)+cp(neigh(k,1,-1,Nc+1)+1)-2*cp(k+1)) ...
            -(lambda*Dt0/2)*cp(k+1)*(cp(k+1)^2-v^2) - epsilon*Dt0/2/v)/dtau;
    end

    %C7. run loop
    for j=1:Nc
        for k=0:(N-1)
            l = j+k*(Nc+1);
            vel(l+1) = vel(l) + dtau*acc(l);
            cp(l+1) = cp(l) + dtau*vel(l+1);%
        end
        for k=0:(N-1)
            l = j+k*(Nc+1);
            acc(l+1) = (1/a^2)*(cp(neigh(l,1,1,Nc+1)+1)+cp(neigh(l,1,-1,Nc+1)+1)-2*cp(l+1)) ...
            -(lambda/2)*cp(l+1)*(cp(l+1)^2-v^2) - epsilon/2/v;    
        end
    end

    %12. combine phi with ap and cp and save combination to file
    tCp = complex(zeros(Tdim,1));
    for j=0:(Tdim-1)
        t = intCoord(j,0,NT);
        x = intCoord(j,1,NT);
        if t<Na
            t = Na - t;
            tCp(j+1) = ap(t+x*(Na+1)+1);
        elseif t<(Na+Nb)
            t = t - Na;
            tCp(j+1) = p(2*(t+x*Nb)+1) + 1i*p(2*(t+x*Nb)+2);
        else
            t = t - Na - Nb;
            tCp(j+1) = cp(t+x*(Nc+1)+1);
        end
    end
    
    tp = vecReal(tCp,Tdim);
    tp(end+1) = p(2*Bdim+1);
    
    
    if 1==1 %loop==0 %printing to terminal
        fprintf('%8s','time', 'runs','N','Nb','L','Lb','R','mass','lambda');
        fprintf('%14s','epsilon','re(action)','im(action)'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('\n');
    end
    fprintf('%8g',toc,runsCount,N,Nb,L,Lb,R,mass,lambda);
    fprintf('%14g%14g%14g\n',epsilon,real(action),imag(action));
    
    actionOut = fopen('data/picAction.dat','a'); %saving action etc to file
    fprintf(actionOut,'%14g',toc,runsCount,d,N,Nb,Na,X,Lb,real(action));
    fprintf(actionOut,'%14g\n',imag(action));
    fclose(actionOut);
    
    save( ['data/picOut',num2str(loop),'.mat'], 'p', 'tCp', 'tp', 'Cp', 'minusDS','DDS','action', 'd', 'N', 'Na', 'Nb' , 'Nc', 'lambda', 'mass', 'R','L','La','Lb','Lc','inP');%saving phi, DDS and minusDS to file
    
end%closing parameter loop

data = load(['data/picOut',num2str(loop),'.mat']);
plotTphi(data);