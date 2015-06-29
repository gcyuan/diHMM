function modelOut = reorderModel(modelIn, oldBinOrderToSorted, oldDomainOrderToSorted)

modelOut = modelIn;

modelOut.emissions = modelIn.emissions(oldBinOrderToSorted,:);
modelOut.transitionB = modelIn.transitionB(oldBinOrderToSorted, oldBinOrderToSorted, oldDomainOrderToSorted);
modelOut.transitionD = modelIn.transitionD(oldDomainOrderToSorted, oldDomainOrderToSorted);
modelOut.initialProbabilities = modelIn.initialProbabilities(oldBinOrderToSorted, oldDomainOrderToSorted);


chrList = getChrList(modelIn.param.genome);
cellTypes = fieldnames(modelIn.states);

for indexCT = 1:length(cellTypes)
    
    for index = 1:length(chrList)
        
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).binStates = ...
            reorderStates(modelIn.states.(cellTypes{indexCT}).(chrList{index}).binStates, oldBinOrderToSorted);
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainStates = ...
            reorderStates(modelIn.states.(cellTypes{indexCT}).(chrList{index}).domainStates, oldDomainOrderToSorted);
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainStatesSmallBins = ...
            reorderStates(modelIn.states.(cellTypes{indexCT}).(chrList{index}).domainStatesSmallBins, oldDomainOrderToSorted);
    end
end

% Add postprocess model

end