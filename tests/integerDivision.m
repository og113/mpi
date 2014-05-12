%program to test the speed of integer divisions

N = uint16(input('input N: '));

tic; %idivide
for j=1:N
    k = idivide(N,j);
end
toc;
disp(['idivide took ',num2str(toc),' seconds']);

tic; %floor
for j=1:N
    k = floor(N/j);
end
toc;
disp(['floor took ',num2str(toc),' seconds']);