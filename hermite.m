hermite%script to leapfrog the euclidean phi solution back through minkowskian
%time - checks on energy conservation are made
%input phi is taken from chosen input file
%using the velocity rather than the position verlet

global d N Nt Ntm NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim;
global R X lambda mass v epsilon theta;

parameters('t');

syms x
roots = vpasolve(x^3 -v^2*x + epsilon/v/lambda,x); %solving V'(p)=0
sort (roots); %roots sorted in ascending order

%1. get ep==ephi from file
%%filename = input('input filename: ','s');
filename = 'data/phiEarly.mat';
data = load(filename);
if isfield(data,'Cp')
    ep = data.Cp;
else
    disp('no vector named Cp in file');
end

%2. initialize p==mphi using last point of ephi and zeros- use complex phi
p = complex(zeros(Mdim));
for j=0:(N-1)
    p(j*Ntm+1) = roots(1); %%ep(j*Nt+1);
end

%3. initialize vel using dS/dp=0 at zero and zeros - complex
vel = complex(zeros(Mdim));
Dt0 = b/2*(-1+1i);
for j=0:(N-1)
    k = j*Ntm;
    vel(k+1) = (Dt0/a^2)*(p(neigh(k,1,1,Ntm)+1)+p(neigh(k,1,-1,Ntm)+1)-2*p(k+1)) ...
        +(lambda*Dt0/2)*p(k+1)*(p(k+1)^2-v^2) + epsilon*Dt0/2/v;
end

%4. initialize acc using phi and expression from equation of motion and zeros-
%complex
acc = complex(zeros(Mdim,1));
dt = -b;
for j=0:(N-1)
    acc(j*Ntm+1) = vel(j*Ntm+1)/dt;
end

%4.5 initialize jerk=dacc/dt
jerk = complex(zeros(Mdim,1));
for j=0:(N-1)
    jerk(j*Ntm+1) = ((Dt0/a^2)*(vel(neigh(k,1,1,Ntm)+1)+vel(neigh(k,1,-1,Ntm)+1)-2*vel(k+1)) ...
        +(lambda*Dt0/2)*vel(k+1)*(3*p(k+1)^2-v^2))/dt;
end

%5. initialize pi, velp, accp and jerkp (where p stands for preliminary) to zeros
velp = complex(zeros(Mdim,1));
accp = complex(zeros(Mdim,1));
jerkp = complex(zeros(Mdim,1));

%6. find expression for E and initialize using phi and zeros - energy
%integrated over time slice, not energy density
H = complex(zeros(Ntm,1));
for j=0:(N-1)
    k = Ntm*j;
    H(1) = H(1) + a*vel(k+1)^2/2 + (p(neigh(k,1,1,Ntm)+1)-p(k+1))^2/a + (a*lambda/8)*(p(k+1)^2-v^2)...
        + a*epsilon*(p(k+1)-v)/2/v;
end

%7. run loop
for j=1:(Ntm-1)
    for k=0:(N-1)
        l=j+k*Ntm;
        p(l) = p(l) + vel(l)*dt + (1/2)*acc(l)*dt^2 + (1/6)*jerk(l)*dt^3; %preliminary step
        velp(l) = vel(l) + acc(l)*dt + (1/2)*jerk(l)*dt^2;
    end
    for k=0:(N-1) %%%%a lot of redundant steps as only need to evaluate neighbours before evaluating acc etc
        l=j+k*Ntm;
        accp(l) = ((Dt0/a^2)*(p(neigh(l,1,1,Ntm)+1)+p(neigh(l,1,-1,Ntm)+1)-2*p(l+1)) ...
        +(lambda*Dt0/2)*pp(l+1)*(p(l+1)^2-v^2) + epsilon*Dt0/2/v)/dt;
        jerkp(l) = ((Dt0/a^2)*(velp(neigh(l,1,1,Ntm)+1)+velp(neigh(l,1,-1,Ntm)+1)-2*velp(l+1)) ...
        +(lambda*Dt0/2)*velp(l+1)*(3*p(l+1)^2-v^2))/dt;
    end
    for k=0:(N-1)
        l=j+k*Ntm;
        vel(l+1) = vel(l) + (1/2)*(acc(l) + accp(l))*dt - (1/12)*(jerkp(l)-jerk(l))*dt^2;
        p(l+1) = p(l) + (1/2)*(vel(l) + vel(l+1))*dt - (1/12)*(accp(l) - acc(l))*dt^2;
    end
    for k=0:(N-1)
        l=j+k*Ntm;
        acc(l+1) = ((Dt0/a^2)*(p(neigh(l,1,1,Ntm)+1)+p(neigh(l,1,-1,Ntm)+1)-2*p(l+1)) ...
        +(lambda*Dt0/2)*p(l+1)*(p(l+1)^2-v^2) + epsilon*Dt0/2/v)/dt;
        H(j+1) = H(j+1) + a*vel(l+1)^2/2 + (p(neigh(l,1,1,Ntm)+1)-p(l+1))^2/a...
            + (a*lambda/8)*(p(l+1)^2-v^2) + a*epsilon*(p(l+1)-v)/2/v;
    end
end

%8. print t and x vs phi
x = xVec(Ntm);
t = real(mTVec);
subplot(2,4,1)
plot3(t,x,real(p),'x')
xlabel('real(t)'), ylabel('x'), zlabel('re(phi)')
subplot(2,4,5)
plot3(t,x,imag(p),'x')
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

%11. print t vs E
shortT = 0:dt:(-Ltm);
subplot(2,4,4)
plot(shortT,real(H),'x')
xlabel('real(t)'), ylabel('real(H)')
subplot(2,4,8)
plot(shortT,imag(H),'x')
xlabel('real(t)'), ylabel('imag(H)')

%12. combine phi with ephi and save combination to file
totalPhi = complex(zeros(Tdim));
for j=0:(Tdim-1)
    if j<Mdim
        totalPhi(j+1) = p(Mdim-j);
    else
        totalPhi(j+1) = ep(Edim-j+Mdim);
    end
end

save data/tp.mat totalPhi;