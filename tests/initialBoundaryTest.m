%program to test whether a periodic instanton in fact satisfies the
%required boundary conditions

global d N Nt Ntm NT NtonN NtmonNt L Lt Ltm a b Edim Mdim Tdim; %defining global variables
global R X lambda mass v epsilon angle theta;

fileNo = input('which data/picOut#.mat file to load? (#) '); %loading periodic instanton
data = load(['data/picOut',num2str(fileNo),'.mat']);

testResult1 = zeros(N,1);

p = data.totalPhi;

for j=0:(N-1)
   testResult1(j+1) = p(2*(j*NT+1)+1)-p(2*j*NT+1);
end

plot(testResult1);