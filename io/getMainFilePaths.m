function file = getMainFilePaths(varargin)
%   GETMAINFILEPATHS Build paths for main files.
%   Input is genome build used. Default is hg19
%
%   Author: Eugenio Marco

if nargin == 0
    refGenome = 'hg19';
else
    refGenome = varargin{1};
end

file.mainDir = pwd;
% All paths relative to main
file.chromSizesDir = 'chromSizes';

file.chromSizesFile = fullfile(file.mainDir, file.chromSizesDir,[refGenome '.txt']);

end
