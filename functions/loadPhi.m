%function to load phi from a .dat file or a .mat file
%argument is filename
function Xphi = loadPhi(filename)
   if ~isempty(strfind(filename, '.dat'))
       inputs = load(filename);
       [m,n]=size(inputs);
       if n==6
           Xphi = inputs(:,5) + 1i*inputs(:,6);
       elseif n==3
           Xphi = inputs(:,3);
       else
           disp('loadPhi error - number of columns in file unnexpected');
       end
   else
       inputs = load(filename);
       Xphi = inputs.Cp;
   end 
end