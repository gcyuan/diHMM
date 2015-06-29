function model = getInitialModelFromModelFinal(modelFinal, param)
%   GETINITIALMODELFROMMODELFINAL, extracts the model parameters from a
%   previous model, so that the EM can be continued from it
%
%   Author: Eugenio Marco

model.markNames = modelFinal.markNames;
model.emissions = modelFinal.emissions;
model.transitionB = modelFinal.transitionB;
model.transitionD = modelFinal.transitionD;
model.initialProbabilities = modelFinal.initialProbabilities;

cellTypes = fields(modelFinal.states)';
chrList = param.trainingChrList;

for index = 1:length(cellTypes)
    for index2 = 1:length(chrList)
        model.states.(cellTypes{index}).(chrList{index2}) =  modelFinal.states.(cellTypes{index}).(chrList{index2});
    end
end

model.binSize = modelFinal.binSize;
model.nDSize = modelFinal.nDSize;
model.logEmissions = modelFinal.logEmissions;
model.logTransitionD = modelFinal.logTransitionD;
model.logTransitionB = modelFinal.logTransitionB;
model.logEmissionsBin2dec = modelFinal.logEmissionsBin2dec;
model.emissionsBin2dec = modelFinal.emissionsBin2dec;

end
