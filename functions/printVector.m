%function to print a vector to a file
function printVector(vectorToPrint,filename)
save(filename,'vectorToPrint','-ascii','-double','-tabs');
end