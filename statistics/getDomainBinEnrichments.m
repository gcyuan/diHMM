function domainBinEnrichments = getDomainBinEnrichments(modelFinal)
%   GETDOMAINBINENRICHMENTS Function to calculate for each domain
%   enrichment into nucleosome-level states
%
%   Author: Eugenio Marco

cellTypes = fieldnames(modelFinal.states);

chrList = fieldnames(modelFinal.states.(cellTypes{1}));

binSize = modelFinal.param.binSize;
nDSize = modelFinal.param.nDSize;

nCT = length(cellTypes);
nB = modelFinal.param.nB;
nD = modelFinal.param.nD;

sizeBinLevelStates = cell(nCT,1);
sizeDomainLevelStates = cell(nCT,1);
domainLevelComposition = cell(nCT,1);
genomeSize= cell(nCT,1);

for indexCT = 1:nCT
    genomeSize{indexCT} = 0; % Should we the same, but just in case
    sizeBinLevelStates{indexCT} = zeros(nB,1);
    sizeDomainLevelStates{indexCT} = zeros(nD,1);
    domainLevelComposition{indexCT} = zeros(nD,nB);
    
    for index =1:length(chrList)
        genomeSize{indexCT} = genomeSize{indexCT} + binSize*length(modelFinal.states.(cellTypes{indexCT}).(chrList{index}).binStates);
        for indexB = 1:nB
            sizeBinLevelStates{indexCT}(indexB) = sizeBinLevelStates{indexCT}(indexB) + ...
                binSize*sum(modelFinal.states.(cellTypes{indexCT}).(chrList{index}).binStates==indexB);
        end
        
        for indexD = 1:nD
            domainLocations = modelFinal.states.(cellTypes{indexCT}).(chrList{index}).domainStates==indexD;
            domainLocationsInBins = reshape(repmat(domainLocations',nDSize,1),...
                length(domainLocations)*nDSize,1);
            sizeDomainLevelStates{indexCT}(indexD) = sizeDomainLevelStates{indexCT}(indexD) + ...
                binSize*nDSize*sum(domainLocations);
            for indexB = 1:nB
                domainLevelComposition{indexCT}(indexD, indexB) = domainLevelComposition{indexCT}(indexD, indexB) +...
                    binSize*sum(modelFinal.states.(cellTypes{indexCT}).(chrList{index}).binStates(domainLocationsInBins)==indexB);
            end
        end
    end
end


cellTypeDomainBinEnrichments = cell(nCT,1);

for indexCT = 1:nCT
    cellTypeDomainBinEnrichments{indexCT} = ...
        domainLevelComposition{indexCT}./repmat(sizeBinLevelStates{indexCT}',modelFinal.param.nD,1)...
        ./repmat(sizeDomainLevelStates{indexCT}/genomeSize{indexCT},1,modelFinal.param.nB);
end

domainBinEnrichments.cellTypeDomainBinEnrichments = cellTypeDomainBinEnrichments;

allCombinedDomainLevelComposition = zeros(nD, nB);
allCombinedSizeBinLevelStates = zeros(nB,1);
allCombinedSizeDomainLevelStates = zeros(nD,1);

for indexCT = 1:nCT
    allCombinedDomainLevelComposition = allCombinedDomainLevelComposition+ domainLevelComposition{indexCT};
    allCombinedSizeBinLevelStates = allCombinedSizeBinLevelStates + sizeBinLevelStates{indexCT};
    allCombinedSizeDomainLevelStates = allCombinedSizeDomainLevelStates + sizeDomainLevelStates{indexCT};
end


allCombinedDomainBinEnrichments = ...
    allCombinedDomainLevelComposition./repmat(allCombinedSizeBinLevelStates',modelFinal.param.nD,1)...
    ./repmat(allCombinedSizeDomainLevelStates/sum(allCombinedSizeDomainLevelStates),1,modelFinal.param.nB);

allCombinedDomainBinEnrichmentsPValues = zeros(size(allCombinedDomainBinEnrichments));

for indexB = 1:nB
    allCombinedDomainBinEnrichmentsPValues(:,indexB) = ...
        getFishersPValues(...
        allCombinedDomainLevelComposition(:,indexB), ...
        sum(allCombinedSizeDomainLevelStates), ...
        allCombinedSizeBinLevelStates(indexB),...
        allCombinedSizeDomainLevelStates);
end
   
domainBinEnrichments.allCombinedDomainBinEnrichments = allCombinedDomainBinEnrichments;

domainBinEnrichments.allCombinedDomainBinEnrichmentsPValues = allCombinedDomainBinEnrichmentsPValues;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pValues = getFishersPValues(hitCounts, totalNBins, totalNHits, binCounts)

pValues = zeros(size(hitCounts));

% Hack for the moment
% Matlab run out of memory calculating the hypergeometric cdf
% We reduce all numbers by a factor redFac
redFac = 100;

for index = 1:length(hitCounts)
    pValues(index) = hygecdf(floor(hitCounts(index)/redFac), floor(totalNBins/redFac), floor(totalNHits/redFac), ...
        floor(binCounts(index)/redFac), 'upper');
end

end

