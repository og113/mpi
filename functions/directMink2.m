<<<<<<< HEAD
%script to directly integrate the euclidean phi solution forwards along minkowskian
%time - checks on energy conservation are made
%input phi is taken from chosen input file

global d N Nt Ntm NT NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
global R X lambda mass v epsilon theta;

parameters('p');
=======
%script to directly integrate the euclidean phi solution back through minkowskian
%time - checks on energy conservation are made
%input phi is taken from chosen input file

global d N Na Nb Nc NT L La Lb Lc a b Adim Bdim Cdim Tdim; %defining global variables
global R X lambda mass v epsilon angle theta;

parameters('p');
Na = Ntm;
Nb = Nt;
Nc = Nt*3;
Cdim = N*Nc;
Lb = Lt;
Lc = b*(Nc-1);
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a

syms x
roots = vpasolve(x^3 -v^2*x + epsilon/v/lambda,x); %solving V'(p)=0
sort (roots); %roots sorted in ascending order

%1. get ep==ephi from file
%%filename = input('input filename: ','s');
<<<<<<< HEAD
data = load('data/tp.mat');
if isfield(data,'totalPhi')
    tPhi = data.totalPhi;
else
    disp('no vector named tPhi in file');
end

%2. initialize mp==mphi using last point of ephi and zeros- use complex phi
mp = complex(zeros(Mdim,1));
for j=0:(N-1)
    mp(j*Ntm+1) = tPhi(j*NT+(NT-1)+1);%%roots(1);%%
=======
data = load('data/picOut0.mat');
if isfield(data,'Cp')
    ep = data.tCp;
else
    disp('no vector named Cp in file');
end

%2. initialize mp==mphi using last point of ephi and zeros- use complex phi
mp = complex(zeros(Cdim,1));
for j=0:(N-1)
    mp(j*Nc+1) = ep((j*NT+NT-1)+1);%%roots(1);%%
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
end

%3. initialize vel - defined at half steps, first step being at t=-1/2,
%vel(t+1/2) := (p(t+1)-p(t))/dt
<<<<<<< HEAD
vel = complex(zeros(Mdim,1));
dtau = b;
Dt0 = dtau;%%b/2*(1-1i*up); - this is surely wrong!!
for j=0:(N-1)
    k = j*Ntm;
=======
vel = complex(zeros(Cdim,1));
dtau = b;
Dt0 = dtau;%%b/2*(-1+1i*up); - this is surely wrong!!
for j=0:(N-1)
    k = j*Nc;
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
    vel(k+1) = 0; %%due to boundary condition
end

%4. initialize acc using phi and expression from equation of motion and zeros-
%complex
<<<<<<< HEAD
acc = complex(zeros(Mdim,1));
for j=0:(N-1)
    k = j*Ntm;
    acc(k+1) = ((Dt0/a^2)*(mp(neigh(k,1,1,Ntm)+1)+mp(neigh(k,1,-1,Ntm)+1)-2*mp(k+1)) ...
=======
acc = complex(zeros(Cdim,1));
for j=0:(N-1)
    k = j*Nc;
    acc(k+1) = ((Dt0/a^2)*(mp(neigh(k,1,1,Nc)+1)+mp(neigh(k,1,-1,Nc)+1)-2*mp(k+1)) ...
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
        -(lambda*Dt0/2)*mp(k+1)*(mp(k+1)^2-v^2) - epsilon*Dt0/2/v)/dtau;
end

%6. find expression for E and initialize using phi and zeros - energy
%integrated over time slice, not energy density
<<<<<<< HEAD
H = complex(zeros(Ntm,1));
for j=0:(N-1)
    k = Ntm*j;
    H(1) = H(1) + a*vel(k+1)^2/2 + (mp(neigh(k,1,1,Ntm)+1)-mp(k+1))^2/a...
=======
H = complex(zeros(Nc,1));
for j=0:(N-1)
    k = Nc*j;
    H(1) = H(1) + a*vel(k+1)^2/2 + (mp(neigh(k,1,1,Nc)+1)-mp(k+1))^2/a...
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
        + (a*lambda/8)*(mp(k+1)^2-v^2) + a*epsilon*(mp(k+1)-v)/2/v;
end

%7. run loop
<<<<<<< HEAD
for j=1:(Ntm-1)
    for k=0:(N-1)
        l = j+k*Ntm;
=======
for j=1:(Nc-1)
    for k=0:(N-1)
        l = j+k*Nc;
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
        vel(l+1) = vel(l) + dtau*acc(l);
        mp(l+1) = mp(l) + dtau*vel(l+1);%
    end
    for k=0:(N-1)
<<<<<<< HEAD
        l = j+k*Ntm;
        acc(l+1) = (1/a^2)*(mp(neigh(l,1,1,Ntm)+1)+mp(neigh(l,1,-1,Ntm)+1)-2*mp(l+1)) ...
        -(lambda/2)*mp(l+1)*(mp(l+1)^2-v^2) - epsilon/2/v;    
        H(j+1) = H(j+1) + a*vel(l+1)^2/2 + (mp(neigh(l,1,1,Ntm)+1)-mp(l+1))^2/a...
=======
        l = j+k*Nc;
        acc(l+1) = (1/a^2)*(mp(neigh(l,1,1,Nc)+1)+mp(neigh(l,1,-1,Nc)+1)-2*mp(l+1)) ...
        -(lambda/2)*mp(l+1)*(mp(l+1)^2-v^2) - epsilon/2/v;    
        H(j+1) = H(j+1) + a*vel(l+1)^2/2 + (mp(neigh(l,1,1,Nc)+1)-mp(l+1))^2/a...
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
            + (a*lambda/8)*(mp(l+1)^2-v^2) + a*epsilon*(mp(l+1)-v)/2/v;
    end
end

%8. print t and x vs phi
<<<<<<< HEAD
x = xVec(Ntm,N);
t = real(mTVec(Ntm));
=======
x = xVec(Nc,N);
t = zeros(Cdim,1);
for j=0:(Cdim-1)
    t(j+1) = b*intCoord(j,0,Nc);
end
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
subplot(2,4,1)
plot3(t,x,real(mp),'x')
xlabel('real(t)'), ylabel('x'), zlabel('re(phi)')
subplot(2,4,5)
plot3(t,x,imag(mp),'x')
xlabel('real(t)'), ylabel('x'), zlabel('imag(phi)')

%9. print t and x vs vel
subplot(2,4,2)
plot3(t,x,real(vel),'x')
xlabel('real(t)'), ylabel('x'), zlabel('re(vel)')
subplot(2,4,6)
plot3(t,x,imag(vel),'x')
xlabel('real(t)'), ylabel('x'), zlabel('imag(vel)')

%10. print t and x vs acc
subplot(2,4,3)
plot3(t,x,real(acc),'x')
xlabel('real(t)'), ylabel('x'), zlabel('re(acc)')
subplot(2,4,7)
plot3(t,x,imag(acc),'x')
xlabel('real(t)'), ylabel('x'), zlabel('imag(acc)')

<<<<<<< HEAD
%11. print t vs E
shortT = 0:dtau:Ltm;
=======
%11. print t vs EMdim = Ntm*N;
shortT = 0:dtau:(Lc);
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
subplot(2,4,4)
plot(shortT,real(H),'x')
xlabel('real(t)'), ylabel('real(H)')
subplot(2,4,8)
plot(shortT,imag(H),'x')
xlabel('real(t)'), ylabel('imag(H)')

%12. combine phi with ephi and save combination to file
<<<<<<< HEAD
totalPhi = complex(zeros(Tdim+Mdim,1));
for j=0:(Tdim+Mdim-1)
    t = intCoord(j,0,2*Ntm+Nt);
    x = intCoord(j,1,2*Ntm+Nt);
    if t<NT
        totalPhi(j+1) = tPhi(t+x*NT+1);
    else
        t = t - NT;
        totalPhi(j+1) = mp(t+x*Ntm+1);
    end
end

save data/tp2.mat totalPhi;
=======
totalPhi = complex(zeros(Tdim,1));
for j=0:(Tdim-1)
    t = intCoord(j,0,NT);
    x = intCoord(j,1,NT);
    if t<(Na+Nb)
        totalPhi(j+1) = ep(t+x*(Na+Nb)+1);
    else
        t = t - Na - Nb;
        totalPhi(j+1) = mp(t+x*Nc+1);
    end
end

save data/tp.mat totalPhi;
>>>>>>> 3da375aa42770427fbfd7b2229902f52969fc49a
