function baseName = getResultsFileBaseName(projectName, resultsType, cellType, runNumber, nB, nD)
%   GETRESULTSFILEBASENAME helper function to get basename
%
% Author: Eugenio Marco

baseName  = fullfile(pwd, 'results',projectName, ...
    ['nB' num2str(nB) '_nD'  num2str(nD)], ['run' num2str(runNumber)],resultsType ,...
    cellType, [cellType '_nD' num2str(nD) '_nB' num2str(nB)]);
end