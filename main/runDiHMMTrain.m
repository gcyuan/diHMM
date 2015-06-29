function runDiHMMTrain(cellTypes, genome, nB, nD, projectName, runNumber, varargin)
%   RUNDIHMMTRAIN Main function to a train a model using diHMM
%
%   Inputs are
%
%   CELLTYPES, cell array with the names of the cell types to use for
%   training the model
%    
%   GENOME, genome to use, needed to get chromosome list and sizes
%
%   NB, number of nucleosome-level states in the model
%
%   ND, number of domain-level states in the model
%
%   PROJECTNAME, name of the folder where the results will be saved in the
%   results folder
%
%   RUNNUMBER, a number, typically 1, to save different results for models
%   trained in different conditions, for example different initializations
%
%   varargin, optional argument, a previous model, which is taken to
%   comtinue training by EM
%
%   For convenience, default values for CELLTYPES, GENOME, NB, ND and
%   GENOME can be set at the beginning of the function
%
%   Author: Eugenio Marco


if ~exist('cellTypes', 'var')
    cellTypes = {'GM12878', 'H1hesc', 'K562'};
end
if ~exist('genome', 'var')
    genome = 'hg19';
end
if ~exist('nB', 'var')
    nB = 30;
end
if ~exist('nD', 'var')
    nD = 30;
end
if ~exist('projectName', 'var')
    projectName = 'modelGHKChr17';
end
if ~exist('runNumber', 'var')
    runNumber = 1;
end

% getTestParameters sets many parameters of the model like the maximum
% number of iteration or initialization method
param = getTestParameters(cellTypes,nB,nD,genome,projectName,runNumber);

% param.doAll controls which chromosome we use to train the model
% by default chr17
switch param.doAll
    case true
        chrList = getChrList(genome);
    otherwise
        chrList = param.chrTodo;
end

param.trainingChrList = chrList;

% read data
for index = 1:length(cellTypes)
    for index2 = 1:length(chrList)
        data.(cellTypes{index}).(chrList{index2}) = readData(cellTypes{index}, chrList{index2}, param);
    end
end

param.nMarks = size(data.(cellTypes{1}).(chrList{1}).binData,2);
param.markNames = data.(cellTypes{1}).(chrList{1}).markNames;

if nargin < 7 % We did not give initialModel
    initialModel = getInitialModel(data, param);
    % initialModel = importChromHMM(param.nB,param);
    
    % Add log of matrices and combine possible emissions into a dec number
    initialModel = postProcessModel(initialModel, param);
    
else % We take the initialModel
    initialModel = getInitialModelFromModelFinal(varargin{1},param);
    clear varargin % Free up memory
end


% call to the main function to train the model
modelTrained = diHMMTrain(data, initialModel, param, chrList);

% Now we load the full genome for decoding
chrList = getChrList(genome);

% read data
for index = 1:length(cellTypes)
    for index2 = 1:length(chrList)
        data.(cellTypes{index}).(chrList{index2}) = readData(cellTypes{index}, chrList{index2}, param);
    end
end

for index = 1:length(cellTypes)
    cellType = cellTypes{index};
    for index2 = 1:length(chrList)
    % NOTE if enough memory the for can be changed to a parfor to decode
    % in parallel
    %     parfor index2 = 1:length(chrList)
        chr = chrList{index2};
        modelIntermediate{index,index2} = getModelIntermediate(cellType,chr,modelTrained,param);
    end
end

modelFinal = modelTrained;
clear modelTrained

for index = 1:length(cellTypes)
    modelFinal.states = rmfield(modelFinal.states,cellTypes{index});
    for index2 = 1:length(chrList)
        modelFinal.states.(cellTypes{index}).(chrList{index2}) = ...
            modelIntermediate{index,index2}.states.(cellTypes{index}).(chrList{index2});
        modelFinal.allLogLiks.(cellTypes{index}).(chrList{index2}) = ...
            modelIntermediate{index,index2}.allLogLiks.(cellTypes{index}).(chrList{index2});
        
    end
end


modelFinal.param = param;

% We assign provisional nucleosome-level and domain-level annotations
% defaults to white states 11 for domains, 13 for nucleosomes

modelFinal = assignAnnotations(modelFinal,'bin',13*ones(1,nB));
modelFinal = assignAnnotations(modelFinal,'domain',11*ones(1,nD));

% save the model
runBaseName = getRunBaseName(projectName,cellTypes, runNumber, nB, nD);
checkBasename(runBaseName)
save([runBaseName '.mat'], 'modelFinal')

% save bed files
for index = 1:length(cellTypes)
    bedFileBaseName = getResultsFileBaseName(projectName, 'bedFiles', cellTypes{index}, runNumber, nB, nD);
%     bedFileBaseName = getBedFileBaseName(projectName, cellTypes{index}, runNumber, nB, nD);
    checkBasename(bedFileBaseName)
    saveBedFiles(bedFileBaseName, modelFinal)
end

% We plot model matrices
plotModel(modelFinal)

end

function modelIntermediate = getModelIntermediate(cellTypes,chrList,modelTrained,param)

data.(cellTypes).(chrList) = readData(cellTypes, chrList, param);
modelIntermediate = diHMMPosteriorStates(data, modelTrained, param);

end





