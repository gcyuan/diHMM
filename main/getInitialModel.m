function model = getInitialModel(data, param)
%   GETINITIALMODEL estimates the initial model parameters using different
%   methods
%   optional methods are
%
%   KCENTER
%   KMEANS
%   RANDOM
%   PCA
%
%   kcenter implemented by Luca Pinello, using function from Yuan Yao PKU
%
%   Other methods implementd by Eugenio Marco

switch param.initialModel
    case 'kcenter'
        % Made by Luca Pinello using
        % kcenter
        % Farthest-First Traversal Algorithm as a 2-approximation for
        % kcenter clustering
        %   Yuan Yao PKU, 2011.02.25
        
        % T.F. Gonzalez. Clustering to minimize the maximum intercluster
        %   distance. Theoretical Computer Science, 38:293-306, 1985.
        cellTypes = fieldnames(data);
        chrList = fieldnames(data.(cellTypes{1}));
        
        model.markNames = param.markNames;
        model.emissions = zeros(param.nB, param.nMarks); % INDEX ORDER: BinState, histone mark
        model.transitionB = zeros(param.nB, param.nB, param.nD); % INDEX ORDER: BinState, BinState, Domain
        model.transitionD = zeros(param.nD, param.nD); % INDEX ORDER: Domain, Domain
        initialProbabilities = ones(param.nB,param.nD);
        model.initialProbabilities = initialProbabilities/sum(initialProbabilities(:));
        
        [allBinData, binDataSizes] = packAllData(data);
        domainDataSizes = binDataSizes/param.nDSize;
        
            
        % BinLevel clustering
        allBinStates = uint8(ones(size(allBinData,1),1)); % Preallocate to state B1
        % We only cluster states with marks
        idxsToClusters = sum(allBinData,2)>0;
        [~, ~, subSolution] = kcenter(double(allBinData(idxsToClusters,:)), param.nB-1);
        allBinStates(idxsToClusters)= uint8(subSolution + 1);
        
        % We cluster domains by the counts they have of binLevel states
        elementaryDomains = reshape(allBinStates,...
            param.nDSize, length(allBinStates)/param.nDSize)';
        elementaryDomainsCounts = levelcounts(ordinal(elementaryDomains,[],1:param.nB),2);
        [~, ~, allDomainStates] = kcenter(elementaryDomainsCounts, param.nD);
        allDomainStatesSmallBins = reshape(repmat(allDomainStates',...
            param.nDSize,1),size(allBinStates));
        % unpack into different cellTypes and chromosomes
        model = unpackAllStates(model,cellTypes, chrList, binDataSizes, domainDataSizes,...
            allBinStates, allDomainStates, allDomainStatesSmallBins);
        
    case 'kmeans'
        cellTypes = fieldnames(data);
        chrList = fieldnames(data.(cellTypes{1}));
        
        model.markNames = param.markNames;
        model.emissions = zeros(param.nB, param.nMarks); % INDEX ORDER: BinState, histone mark
        model.transitionB = zeros(param.nB, param.nB, param.nD); % INDEX ORDER: BinState, BinState, Domain
        model.transitionD = zeros(param.nD, param.nD); % INDEX ORDER: Domain, Domain
        initialProbabilities = ones(param.nB,param.nD);
        model.initialProbabilities = initialProbabilities/sum(initialProbabilities(:));
        
        [allBinData, binDataSizes] = packAllData(data);
        domainDataSizes = binDataSizes/param.nDSize;
        allDomainData = downSampleBin(allBinData, param.nDSize);
        allBinStates = uint8(kmeans(double(allBinData), param.nB,...
            'start','cluster','emptyaction','singleton'));
        allDomainStates = uint8(kmeans(allDomainData, param.nD,...
            'start','cluster','emptyaction','singleton'));
        allDomainStatesSmallBins = reshape(repmat( allDomainStates',...
            param.nDSize,1),size(allBinStates));
        model = unpackAllStates(model,cellTypes, chrList, binDataSizes, domainDataSizes,...
            allBinStates, allDomainStates, allDomainStatesSmallBins);
    case 'random'
        error('Broken for the moment')
        % We set random initial parameters, with the option of using a seed
        
        model.emissions = rand(param.nB,param.nMarks); % INDEX ORDER: BinState, histone mark
        model.transitionB = repmat(initialMatrixValues(param.nB), [1 1 param.nD]); % INDEX ORDER: BinState, BinState, Domain
        model.transitionD = initialMatrixValues(param.nD); % INDEX ORDER: Domain, Domain
        initialProbabilities = ones(param.nB,param.nD);
        model.initialProbabilities = initialProbabilities/sum(initialProbabilities(:));
        
        model.transitionD = modifyMatrix(model.transitionD,param.nTimes);
        for index = 1:param.nD
            model.transitionB(:,:,index) = modifyMatrix(model.transitionB(:,:,index),param.nTimes);
        end
        
    case 'PCA'
        chrList = fieldnames(data);
        model.markNames = param.markNames;
        model.emissions = zeros(param.nB, param.nMarks); % INDEX ORDER: BinState, histone mark
        model.transitionB = zeros(param.nB, param.nB, param.nD); % INDEX ORDER: BinState, BinState, Domain
        model.transitionD = zeros(param.nD, param.nD); % INDEX ORDER: Domain, Domain
        initialProbabilities = ones(param.nB,param.nD);
        model.initialProbabilities = initialProbabilities/sum(initialProbabilities(:));
        
        for index =1:length(chrList)
            % uint8, works for param.nD < 255
            %             model.states.(chrList{index}).domainStates = uint8(kmeans(domainData, param.nD,...
            %                 'start','cluster','emptyaction','singleton'));
            % uint8, works for param.nB < 255
            model.states.(chrList{index}).binStates = uint8(kmeans(double(data.(chrList{index}).binData), param.nB,...
                'start','cluster','emptyaction','singleton'));
            elementaryDomains = reshape(model.states.(chrList{index}).binStates,...
                param.nDSize, length(model.states.(chrList{index}).binStates)/param.nDSize)';
            elementaryDomainsCounts = levelcounts(ordinal(elementaryDomains,[],1:param.nB),2);
            
            [~, elementaryDomainsCountsScore] = princomp(elementaryDomainsCounts);
            
            model.states.(chrList{index}).domainStates = uint8(kmeans(elementaryDomainsCountsScore, param.nD,...
                'start','cluster','emptyaction','singleton'));
            
            model.states.(chrList{index}).domainStatesSmallBins = ...
                reshape(repmat( model.states.(chrList{index}).domainStates',...
                param.nDSize,1),size(model.states.(chrList{index}).binStates));
        end
end
model = updateInitialModelParameters(model, data, param);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matrix = initialMatrixValues(nB)
% We split 90% for the diagonals and the rest evenly split, with sums for
% rows equal to 1
defaultDiagonal = 0.9;
offDiagonal = (1-defaultDiagonal)/(nB-1);
matrix = offDiagonal*ones(nB);
dummyDiagonal = diag(repmat(defaultDiagonal-offDiagonal,1,nB));
matrix = matrix + dummyDiagonal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matrix = modifyMatrix(matrix,nTimes)
% We modify randomly the initial matrix by affecting two elements, but
% preserving row sum to 1
nCols = size(matrix,1);
for index = 1:nTimes
    % We choose two columns and a row randomly
    cols = floor(1+nCols*rand(1,2));
    while cols(2)==cols(1)
        cols(2) = floor(1+nCols*rand(1));
    end
    row = floor(1+nCols*rand);
    subMatrix = matrix(row,cols);
    smallestDistance = min(min(1-subMatrix(:)),min(subMatrix(:)));
    step = smallestDistance/20;
    randomSign = sign(randn);
    subMatrix = subMatrix + randomSign* [step -step];
    matrix(row,cols) = subMatrix;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From Patrick Mineault
% http://xcorr.net/2008/03/13/binning-in-matlab-a-one-liner/
% modified to crop end, fixed dimnum to 1
function outv = downSampleBin(invec,amount)
dimnum = 1;
invec = invec(1:end-mod(size(invec,1),amount),:);
sizebefore = size(invec);
sizemiddle = [sizebefore(1:dimnum-1),amount,sizebefore(dimnum)/amount,sizebefore(dimnum+1:length(sizebefore))];

sizeafter = sizebefore;

sizeafter(dimnum) = sizeafter(dimnum)/amount;

outv = reshape(sum(reshape(invec,sizemiddle),dimnum),sizeafter);
end

function [allData, dataSizes] = packAllData(data)
% Function to pack together all the data
cellTypes = fieldnames(data);
chrList = fieldnames(data.(cellTypes{1}));
dataSizes = zeros(1,length(cellTypes)*length(chrList));

nCols = size(data.(cellTypes{1}).(chrList{1}).binData,2);
n = 1;
% Preallocation
for index =1:length(cellTypes)
    for index2 =1:length(chrList)
        dataSizes(n) = size(data.(cellTypes{index}).(chrList{index2}).binData,1);
        n = n+1;
    end
end
nRows = sum(dataSizes);
allData = false(nRows,nCols);
n = 1;
limits = [1 cumsum(dataSizes)+1];
for index =1:length(cellTypes)
    for index2 =1:length(chrList)
        allData(limits(n):(limits(n+1)-1),:) = data.(cellTypes{index}).(chrList{index2}).binData;
        n = n+1;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = unpackAllStates(model,cellTypes, chrList, binDataSizes, domainDataSizes,...
    allBinStates, allDomainStates, allDomainStatesSmallBins)
% Function to unpack all the states into cellTypes and chromosomes

limitsBins = [1 cumsum(binDataSizes)+1];
limitsDomains = [1 cumsum(domainDataSizes)+1];
n=1;

for index =1:length(cellTypes)
    for index2 =1:length(chrList)
        % uint8, works for param.nD < 255
        model.states.(cellTypes{index}).(chrList{index2}).domainStates = ...
            allDomainStates(limitsDomains(n):(limitsDomains(n+1)-1));
        % uint8, works for param.nB < 255
        model.states.(cellTypes{index}).(chrList{index2}).binStates = ...
            allBinStates(limitsBins(n):(limitsBins(n+1)-1));
        model.states.(cellTypes{index}).(chrList{index2}).domainStatesSmallBins = ...
            allDomainStatesSmallBins(limitsBins(n):(limitsBins(n+1)-1));
        n = n + 1;
    end
end


end
