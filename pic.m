%script to find the periodic instanton of the stable phi^4 with negative
%mass
clear all;

global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim; %defining global variables
global R epsilon dE minima angle amp A;

date = '17.9.14';

aq.inputP = 'b'; %struct to hold answers to questions aq short for 'answers to questions' - defauLbs in initialization
aq.pot = 1;
aq.perturbResponse = 'n';
aq.loopResponse = 'n';
aq.parameterChoice = 'N';
aq.minValue = 32;
aq.maxValue = 64;
aq.totalLoops = 5;

inP = aq.inputP; %just because I write this a lot
parameters(inP,aq.pot); %assigning global variables according to parameters.m

fprintf('%10s','N', 'Na','Nb','Nc','L','Lb','R','epsilon','dE'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
fprintf('\n');
fprintf('%10g',N,Na,Nb,Nc,L,Lb,R,epsilon);
fprintf('%10g\n',dE);

tempAq = askQuestions; %asking questions
fields = fieldnames(tempAq);
for j=1:numel(fields) %using non-empty responses to questions to fill in aq
    if ~isempty(tempAq.(fields{j}))
       aq.(fields{j}) = tempAq.(fields{j}); 
    end
end
inP = aq.inputP; %just because I write this a lot
parameters(inP,aq.pot);

negVal = 0;
negVec = zeros(2*N*Nb+1,1);
if strcmp(inP,'q') || strcmp(inP,'p') || strcmp(inP,'i') || 1==1
    eigenData = load(['data/',date,'/eigens.mat']);
    negVal = eigenData.D;
    negVec = eigenData.V;
end
if strcmp(inP,'i')
    inputData = load(['data/',date,'/picOut0.mat']);
end

clear x epsi;
%assigning choice of potential
if aq.pot==1
    V = @(x) (1.0/8.0)*(x.^2-1.0).^2 - (epsilon/2.0)*(x-1.0);
    dV = @(x) (x.*(x.^2 - 1.0))/2 - epsilon/2.0;
    ddV = @(x) (3*x.^2 - 1.0)/2;
    toEdge3 = 3e-1;
    toEdge1 = 3e-1;
elseif aq.pot==2
    W = @(x) exp(-x.^2).*(x + x.^3 + x.^5);
    dW = @(x) exp(-x.^2).*(- 2*x.^6 + 3*x.^4 + x.^2 + 1);
    ddW = @(x) 2*x.^3*exp(-x.^2).*(2*x.^4 - 9*x.^2 + 5);
    V = @(x) (1/2)*(x+1).^2.*(1-epsilon*W((x-1)/A));
    dV = @(x) (x+1).*(1-epsilon.*W((x-1)/A)) - (1/2)*(x+1).^2.*(epsilon/A)*dW((x-1)/A);
    ddV = @(x) (1-epsilon*W((x-1)/A)) - (x+1).*(epsilon/A)*dW((x-1)/A) + (1/2)*(x+1).^2*(epsilon/A^2).*ddW((x-1)/A);
    if strcmp(inP,'b')
        integrand = @(x) (2.0*V(x)).^(-0.5);
        toEdge3 = 5e-1; tempToEdge3 = toEdge3;
        toEdge1 = 1e-1;
        minRho = integral(integrand,minima(3)-toEdge3,0);
        for y=tempToEdge3:(-1e-2):1e-2
            smallerRho = integral(integrand,minima(3)-y,0);
            if abs(imag(smallerRho))<eps
                minRho = smallerRho;
                toEdge3 = y;
            end
        end
        maxRho = integral(integrand,minima(1)+toEdge1,0);
        if abs(imag(maxRho))>eps || abs(imag(minRho))>eps
            disp('maxRho or minRho is complex, consider changing integration points');
            return
        else
            phiRho = (minima(3))-toEdge3:-1e-2:(minima(1)+toEdge1);
            intHand = @(x) integral(integrand,x,0);
            Rho = arrayfun(intHand,phiRho);
        end
    end
else
    disp('choice of potential not available');
end

for loop=0:(aq.totalLoops-1) %starting parameter loop, note: answ.totalLoops=1 if answ.loopResponse='n'
    if strcmp(aq.loopResponse,'y') && strcmp(aq.parameterChoice,'N')
        N = aq.minValue + floor(loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1));
        changeParameters(N,'N',inP,aq.pot);
    elseif strcmp(aq.loopResponse,'y')
        loopParameter = aq.minValue + loop*(aq.maxValue - aq.minValue)/(aq.totalLoops-1);
        changeParameters (loopParameter,aq.parameterChoice, inP,aq.pot);
    end
    if strcmp(inP,'i') && loop>0
        inputData = load(['data/',date,'/picOut',num2str(loop-1),'.mat']);
    end
    tic; %starting the clock
    
    %S1 = 2/3; %this is twice the value in the coleman paper, in the thin wall limit
    clear x;
    integrandS1 = @(x) (2.0*V(x)).^0.5;
    S1 = integral(integrandS1,minima(1)+toEdge1,minima(3)-toEdge3);
    if aq.pot==1
        S1 = 2/3;
    end
    if abs(imag(S1))>eps
        disp('S1 integral went too clse to boundaries');
    end
    twAction = -solidAngle(d)*dE*R^d/d + solidAngle(d)*R^(d-1)*S1; %thin-wall bubble action
    M = 2.0/3.0; %soliton mass, with m^2/lambda factored off and in units of m
    betaL = R/3; %determines range over which tanh(x) is used
    betaR = R/3;
    if aq.pot==2
        fitEnd = max(floor(length(Rho)/24),6);
        clear x;
        fitobjectL = fit(Rho(1:fitEnd)',phiRho(1:fitEnd)','poly1');
        betaL = min(abs((minima(3)-phiRho(1))/fitobjectL.p1),5*a);
        fitobjectR = fit(Rho(end-fitEnd:end)',phiRho(end-fitEnd:end)','poly1');
        betaR = min(abs((minima(1)-phiRho(end))/fitobjectR.p1),5*a);
    end
    if ~strcmp(aq.parameterChoice,'amp')
        amp = 2*(Lb-R)/R; %determines admixture of negative mode - trial and error
    end
    action = complex(2);
    
    actionLast = complex(twAction); %defining some quantities to stop Newton-Raphson loop when action stops varying
    runsCount = 0;
    actionTest = 1;
    deltaTest = 1;
    normTest = 1;
    maxTest = 1;
    gaussianTest = 1;
    closenessA = 1;
    closenessD = 1;
    closenessN = 1e-5;
    closenessM = 1e-4;
    minRuns = 3;
    
    p = zeros(2*Bdim+1,1); %phi, in the euclidean domain
    
    for j=0:(Bdim-1) %assigning input phi according to inputP, note that the periodic instanton input is explicitly 2d
        t = eCoord(j,0);
        x = eCoord(j,1);
        clear arg;
        rhoH = @(arg1,arg2) (-arg1.^2 + arg2.^2).^0.5;
        rho = real(rhoH(t,x)); %rho should be real even without the real()
        rho1 = real(rhoH(t,x+R*cos(angle)));
        rho2 = real(rhoH(t,x-R*cos(angle))); 
        if R<betaL || R<betaR
            disp(['R is too small. not possible to give thinwall input. it should be more that ',num2str(max(betaL,betaR))]);
            break
        else
            if strcmp(inP,'t')
                p(2*j+1) = minima(3);
            elseif strcmp(inP,'f')
                p(2*j+1) = minima(1);
            elseif strcmp(inP,'b') || strcmp(inP,'q')
                if aq.pot==1
                    if rho<(R-betaL)
                        p(2*j+1) = minima(3);
                    elseif rho>(R+betaR)
                        p(2*j+1) = minima(1);
                    else
                        p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(1)-minima(3))*tanh((rho-R)/2)/2.0;
                    end
                elseif aq.pot==2
                    if (rho-R)>=minRho && (rho-R)<=maxRho
                        [temp,rhoPos] = min(abs((rho-R)-Rho));
                        p(2*j+1) = phiRho(rhoPos);
                    elseif (rho-R)>=(minRho-betaL) && (rho-R)<maxRho %just to smooth edges
                        p(2*j+1) = fitobjectL(rho-R);
                    elseif (rho-R)>minRho && (rho-R)<=(maxRho+betaR) %just to smooth edges
                        p(2*j+1) = fitobjectR(rho-R);
                    elseif (rho-R)<(minRho-betaL)
                        p(2*j+1) = minima(3);    
                    else
                        p(2*j+1) = minima(1);
                    end
                end
            elseif strcmp(inP,'p')
                if rho1<(R-betaL) && rho2<(R-betaR)
                    p(2*j+1) = minima(3);
                elseif rho1>(R+betaL) || rho2>(R+betaR)
                    p(2*j+1) = minima(1);
                elseif real(eCoord(j,1))>0 %explicitly 2d here, note that the coord should be real
                    p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(1)-minima(3))*tanh((rho1-R)/2)/2.0;
                    %%pNeg(2*j+1) = v/cosh(mass*(rho1-R)/2)^2;
                elseif real(eCoord(j,1))<0
                    p(2*j+1) = (minima(1)+minima(3))/2.0 + (minima(1)-minima(3))*tanh((rho2-R)/2)/2.0;
                    %%pNeg(2*j+1) = v/cosh(mass*(rho2-R)/2)^2;
                else
                    p(2*j+1) = (minima(1)+minima(3))/2.0; %if eCoord(j,1) == 0
                end
            elseif strcmp(inP,'i')
                p(2*j+1) = inputData.p(2*j+1);
            end
            if strcmp(inP,'q') || ((strcmp(inP,'i') && strcmp(inputData.inP,'b')))
                p(2*j+1) = p(2*j+1) + amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+1); %adding the negative eigenvector
                p(2*j+2) = p(2*j+2) + amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+2);
            end
        end
    end
    
    p(2*Bdim+1) = 0.5; %initializing Lagrange parameter for dp/dx zero mode
    
    if (strcmp(inP,'p') || strcmp(inP,'q')) && 1 %fixing input periodic instanton to have zero time derivative at time boundaries
        open = 1.0; %value of 0 assigns all weight tpico boundary, value of 1 to neighbour of boundary
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
        if norm(Chi0) < eps
            disp('norm(Chi0)<eps, done nothing about it');
        end
        
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
        pot = complex(0);
        clear c3;
        
        for j=0:(Bdim-1)
            t = intCoord(j,0,Nb);
            x = intCoord(j,1,Nb);
            dtj = dt(j);
            Dtj = Dt1(j);
            siteMeasure = a*Dtj; %for sites in time
            linkMeasure = a*dtj; %for links in time
            
            if abs(Chi0(j+1))>eps %zero mode lagrange constraint
                DDSm(c3) = 2*j+1; DDSn(c3) = 2*Bdim+1; DDSv(c3) = a*Chi0(j+1);
                DDSm(c3) = 2*Bdim+1; DDSn(c3) = 2*j+1; DDSv(c3) = a*Chi0(j+1);  
                minusDS(2*j+1) = minusDS(2*j+1) - a*Chi0(j+1)*p(2*Bdim+1);
                minusDS(2*Bdim+1) = minusDS(2*Bdim+1) - a*Chi0(j+1)*p(2*j+1);                      
            end
            
            pot = pot + siteMeasure*V(Cp(j+1));
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
                    temp1 = siteMeasure*(2*Cp(j+1)/a^2 + dV(Cp(j+1)));
                    temp2 = siteMeasure*(2/a^2 + ddV(Cp(j+1)));
                    
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
        action = kinetic - pot;
        DDSm(DDSv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
        DDSn(DDSv==0) = [];
        DDSv(DDSv==0) = [];
        DDS = sparse(DDSm,DDSn,DDSv,2*Bdim+1,2*Bdim+1);
        
        small = norm(minusDS); %normalising problem
        smaller = small/(2*Bdim+1);
        cutoff = eps;
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
        if undetermined>eps
            fprintf('%12s','size - structural rank = ');
            fprintf('%12g\n',undetermined);
            break %goes to beginning of while loop and starts again
        end
        
        limit = 1e15*min([closenessA,closenessD,closenessN,closenessM]); %if c>limit numerical errors will be significant
        c = condest(DDS); %finding estimate and bound for condition number
        if c>limit
            fprintf('%12s','condest = ');
            fprintf('%12g\n',c);
        end
        
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
        diff = DDS*delta-minusDS;
        gaussianTest = [gaussianTest, max(abs(diff))];
        
        saveEarly = ['data/',date,'/picEarly',inP,num2str(loop),num2str(runsCount),'.mat'];
        save(saveEarly);
        
        if size(delta,1)==1
            p = p + delta'; %p -> p'
        else
            p = p + delta;
        end

        %pNeg and pZer plus log(det(DDS)) stuff

        stopTime = toc;
        [Xwait, Xaq] = convergenceQuestions(aq, runsCount, stopTime, action,gaussianTest); %discovering whether or not n-r has converged, and stopping if it is wildly out
        aq = Xaq;
        
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
            -Dt0*dV(ap(k+1)))/dtau;
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
            -dV(ap(l+1));
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
            -Dt0*dV(cp(k+1)))/dtau;
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
            -dV(cp(l+1));    
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
            t = t - Na - Nb + 1;
            tCp(j+1) = cp(t+x*(Nc+1)+1);
        end
    end
    
    tp = vecReal(tCp,Tdim);
    tp(end+1) = p(2*Bdim+1);
    
    
    if 1==1 %loop==0 %printing to terminal
        fprintf('%8s','time', 'runs','N','Nb','L','Lb','R');
        fprintf('%14s','dE','re(action)','im(action)'); %can add log|det(DDS)| and 0-mode and neg-mode etc.
        fprintf('\n');
    end
    fprintf('%8g',toc,runsCount,N,Nb,L,Lb,R);
    fprintf('%14g%14g%14g\n',dE,real(action),imag(action));
    
    actionOut = fopen(['data/',date,'/picAction',inP,'.dat'],'a'); %saving action etc to file
    fprintf(actionOut,'%14g',toc,runsCount,d,N,Nb,Na,R,epsilon,Lb,real(action));
    fprintf(actionOut,'%14g\n',imag(action));
    fclose(actionOut);
    
    saveEnd = ['data/',date,'/picOut',inP,num2str(loop),'.mat'];
    save(saveEnd);
    
    if strcmp(inP,'b')
        [V,D] = eigs(DDS,1,-10);
        save (['data/',date,'/eigens.mat'],'D','V');
    end
    
end%closing parameter loop

data = load(saveEnd);
plotTphi(data);