%script to carry out a leap frog integration of second order equation of
%motion - using the harmonic oscillator
%using the velocity (rather than the position) verlet

L = 20; %length of time
N = 100; %number of points
dt = L/(N-1); %time step
t = 0:dt:L;

x = zeros(N); %initializing vectors
v = zeros(N);
vh = zeros(N); %velocity at half time steps
a = zeros(N);
H = zeros(N); %energy

x(1) = input('initial position: '); %initial condition
v(1) = input('initial velocity: ');
a(1) = -x(1);
H(1) = v(1)^2/2 + x(1)^2/2;

for j=1:(N-1)
    vh(j) = v(j) + (dt/2)*a(j);
    x(j+1) = x(j) + dt*vh(j);
    a(j+1) = -x(j+1);
    v(j+1) = vh(j) + (dt/2)*a(j+1);
    H(j+1) = v(j+1)^2/2 + x(j+1)^2/2;
end

subplot(2,2,1)
plot(t,x)
xlabel('t'),ylabel('x')
subplot(2,2,2)
plot(t,v)
xlabel('t'), ylabel('v')
subplot(2,2,3)
plot(x,v)
xlabel('x'),ylabel('v')
subplot(2,2,4)
plot(t,H)
xlabel('t'), ylabel('H')
