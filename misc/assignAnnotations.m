function modelOut = assignAnnotations(modelIn, annotationType, assignedAnnotations)
% ASSIGNANNOTATIONS helper function to assign annotations to bins and
% domains
%
%   Author: Eugenio Marco

switch annotationType
    case 'bin'
        modelIn.assignedBinAnnotations = assignedAnnotations;
    case 'domain'
        modelIn.assignedDomainAnnotations = assignedAnnotations;       
end

modelOut = modelIn;

end