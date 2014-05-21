%constructs relevant matrix and checks largest modulus eigenvalues, which
%are relevant to stability

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
p = complex(zeros(Mdim,1));
for j=0:(N-1)
    p(j*Ntm+1) = roots(1);%%ep(j*Nt+1);
end

%3 initialize vectors of sparse matrix elements
nonzeros = 3*N;
Mm = zeros(nonzeros,1);
Mn = zeros(nonzeros,1);
Mv = complex(zeros(nonzeros,1));

%4. stating some parameters
bonc = 1;
dt = -b/bonc;
Dt0 = dt;%%%b/2*(-1/bonc+1i);

clear c3;

for j=0:(N-1)
    Mm(c3) = j+1; Mn(c3) = j+1; Mv(c3) = 1 + dt*(-2*(Dt0/a^2)...
        +(Dt0*lambda/2)*(3*p(j*Ntm+1)^2 - v^2));
    if j==0
        Mm(c3) = j+1; Mn(c3) = j+2; Mv(c3) = dt*Dt0/a^2;
        Mm(c3) = j+1; Mn(c3) = N; Mv(c3) = dt*Dt0/a^2;
    elseif j==(N-1)
        Mm(c3) = j+1; Mn(c3) = 1; Mv(c3) = dt*Dt0/a^2;
        Mm(c3) = j+1; Mn(c3) = j; Mv(c3) = dt*Dt0/a^2;
    else
        Mm(c3) = j+1; Mn(c3) = j+2; Mv(c3) = dt*Dt0/a^2;
        Mm(c3) = j+1; Mn(c3) = j; Mv(c3) = dt*Dt0/a^2;
    end
end

clear c3;

Mm(Mv==0) = []; %dropping unused nonzero elements of DDS (i.e. zeros)
Mn(Mv==0) = [];
Mv(Mv==0) = [];
M = sparse(Mm,Mn,Mv,N,N);

d = abs(eigs(M,1))