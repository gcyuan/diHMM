function changeAnnotations(cellTypes, nB, nD, projectName, runNumber,...
    annotationType, stateNumbers, newAnnotations)
%   CHANGEANNOTATIONS Function to change nucleosome-level or domain-level
%   annotations 
%
%   ANNOTATIONTYPE can be bin (for nucleosome) or domain
%
%   STATENUMBERS are the states we want to change their annotation
%
%   NEWANNOTATIONS are the new nucleosome or domain annotations. To see
%   available annotations run PLOTDOMAINANNOTATIONS, or PLOTNUCLEOSOMEANNOTATIONS
%
%   Author: Eugenio Marco

baseName = getRunBaseName(projectName, cellTypes, runNumber, nB, nD);
load([baseName '.mat'])

% Extract current order
[~, newBinOrder] = sort(modelFinal.assignedBinAnnotations);
[~, newDomainOrder] = sort(modelFinal.assignedDomainAnnotations);

switch annotationType
    case 'bin'
        modelFinal.assignedBinAnnotations(newBinOrder(stateNumbers)) = newAnnotations;
    case 'domain'
        modelFinal.assignedDomainAnnotations(newDomainOrder(stateNumbers)) = newAnnotations;       
end

save([baseName '.mat'], 'modelFinal')

end