function param = getTestParameters(cellTypes,nB,nD,genome,projectName,...
    runNumber,binSize)
%   GETTESTPARAMETERS sets many of the parametes of the model
%   by default BINSIZE is 200 bases
%
%   DOALL controls whether to train the model using all chromosomes,
%   default false
%
%   CHRTODO controls the chromosome to train the model when DOALL is false
%   default is chr17
%
%   Author: Eugenio Marco

if ~exist('binSize', 'var')
    binSize = 200;
end

% default parameters for the moment
param.dataLocation = 'data'; % relative to main folder
param.cellTypes = cellTypes;
% param.cellTypes = 'H1hesc';
param.binSize = binSize;
param.binDataFolder = ['binarized' num2str(param.binSize) 'bpFiles'];
param.bin2decDataFolder = ['bin2dec' num2str(param.binSize) 'bpFiles'];

% Define filenames
param.genome = genome;
param.projectName = projectName;
param.runNumber = runNumber;

param.filesMain = getMainFilePaths(param.genome); % Default is now hg19, use argument for another one

param.chromSizes = getChromSizes(param.filesMain.chromSizesFile);
param.maxIter = 500; % max number of EM iterations
param.tol = 1e-6; % tolerance

param.doAll = false;
param.chrTodo = {'chr17'}; % Contains Hoxb1 in human

for index = 1:length(cellTypes)
    param.fullBinDataFolder.(cellTypes{index}) =  fullfile(param.filesMain.mainDir,param.dataLocation,...
        param.genome, cellTypes{index},param.binDataFolder);
    param.fullBin2decDataFolder.(cellTypes{index}) =  fullfile(param.filesMain.mainDir,param.dataLocation,...
        param.genome, cellTypes{index},param.bin2decDataFolder);
end


param.nD = nD; % Number of domains
param.nDSize = 20; % Minimum domain size, as number of bins.
param.nB = nB; % Number of binLevel states

% param.initialModel = 'random';
% param.initialModel = 'kmeans';
% param.initialModel = 'hierarchical';
param.initialModel = 'kcenter';


param.useRandomizeInitialSeed = 1;
param.randomizeInitialSeed = sum(100*clock);
param.initialSeed = 12;
param.nTimes = 5; % times we randomize each initial conditions matrices

if param.useRandomizeInitialSeed
    stream0 = RandStream('mt19937ar','Seed',param.randomizeInitialSeed);
    RandStream.setGlobalStream(stream0);
else
    % We choose the same seed for debugging purposes
    stream0 = RandStream('mt19937ar','Seed',param.initialSeed);
    RandStream.setGlobalStream(stream0);
end

param.emissionsKick = .01; % Used to remove zeros in initial estimates

param.calculatePosteriorAndMarginal = false;
param.saveBedgraphPosterior = false;

