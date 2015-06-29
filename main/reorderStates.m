function statesOut = reorderStates(statesIn, oldOrderToSorted)
%   REORDERSTATES Function to reorder the states
%
%   Author: Eugenio Marco

statesOut = statesIn;

n = length(oldOrderToSorted);

for index = 1:n
    statesOut(statesIn ==  oldOrderToSorted(index)) = index;
    
end

end