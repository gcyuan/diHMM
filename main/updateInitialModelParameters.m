function modelOut = updateInitialModelParameters(modelIn, data, param)
%   UPDATEINITIALMODELPARAMETERS
%   Given a state sequence, calculate the emission and transition
%   probabilities
%
%   Author: Eugenio Marco

cellTypes = fieldnames(data);
chrList = fieldnames(data.(cellTypes{1}));

modelOut = modelIn;

% Domain-level transitions
for index = 1:length(cellTypes)
    for index2 = 1:length(chrList)
        for indexD1 = 1:param.nD
            for indexD2 = 1:param.nD
                modelOut.transitionD(indexD1,indexD2) = modelOut.transitionD(indexD1,indexD2) +...
                    getNTransitions(modelIn.states.(cellTypes{index}).(chrList{index2}).domainStates,indexD1,indexD2);
            end
        end
    end
end

% Normalization, we add eps to the denominator
modelOut.transitionD(modelOut.transitionD==0) = min(modelOut.transitionD(modelOut.transitionD~=0));

modelOut.transitionD = modelOut.transitionD./(repmat(sum(modelOut.transitionD,2),1,param.nD)+eps);

% Bin-level transitions
for index = 1:length(cellTypes)
    for index2 = 1:length(chrList)
        for indexD = 1:param.nD
            selectedDomain = modelIn.states.(cellTypes{index}).(chrList{index2}).domainStatesSmallBins == indexD;
            % Given a domain, we select the bin states directly underneath +
            % the state right before the domain (circshift part)
            % We send that to getNTransitions for all possible binState
            % combinations
            stateList = double(modelIn.states.(cellTypes{index}).(chrList{index2}).binStates(1:end-1)).*...
                (selectedDomain(1:end-1)|circshift(selectedDomain(1:end-1),-1));
            for indexB1 = 1:param.nB
                for indexB2 = 1:param.nB
                    modelOut.transitionB(indexB1, indexB2, indexD) = ...
                        modelOut.transitionB(indexB1, indexB2, indexD) +...
                        getNTransitions(stateList,indexB1,indexB2);
                end
            end
        end
    end
end

% Normalization, we add eps to the denominator
% We also remove zeros by setting them equal to the minimum value
modelOut.transitionB(modelOut.transitionB==0) = min(modelOut.transitionB(modelOut.transitionB~=0));

for indexD = 1:param.nD
    modelOut.transitionB(:,:,indexD) = modelOut.transitionB(:,:,indexD)./...
        (repmat(sum(modelOut.transitionB(:,:,indexD),2),1,param.nB)+eps);
end


% Emissions
for indexB =1:param.nB
    countNBinStates = 0;
    countEmissionsBinStates = 0;
    for index = 1:length(cellTypes)
        for index2 = 1:length(chrList)
            states = modelIn.states.(cellTypes{index}).(chrList{index2}).binStates==indexB;
            countNBinStates = countNBinStates + ...
                sum(states);
            countEmissionsBinStates = countEmissionsBinStates + ...
                sum(data.(cellTypes{index}).(chrList{index2}).binData(states,:),1);
            
            %             modelOut.emissions(indexB,:) = ...
            %                 mean(data.(chrList{index}).binData(modelIn.states.(chrList{index}).binStates==indexB,:));
        end
    end
    modelOut.emissions(indexB,:) = countEmissionsBinStates / countNBinStates;
end

modelOut.emissions = param.emissionsKick + (1 - param.emissionsKick) * modelOut.emissions;

modelOut.binSize = param.binSize;
modelOut.nDSize = param.nDSize;


% logOneMinusEmissions = log(1-modelOut.emissions+eps); % + eps to avoid Inf, probably unnecesary
% logEmissions = log(modelOut.emissions+eps);

% We prefill the combination matrix with NaNs and fill it as necessary
% For a given dimension element 1 is emission 0, element 2 is emission 1
% combinationsBinStateEmissions = nan([repmat(2,1,param.nMarks) param.nB]);
% for index = 1:length(chrList)
%     nBins = size(data.(chrList{index}).binData,1);
%     logEmissionsBin.(chrList{index}) = zeros(param.nB,nBins);
%     for index2 = 1:nBins
%         indsObs = num2cell(data.(chrList{index}).binData(index2,:)+1);
%         indsObsFirstState = [indsObs {1}];
%         indsObsAllStates = [indsObs {1:param.nB}];
%         
%         %     value = squeeze(combinationsBinEmissions(inds{:}))
%         if isnan(combinationsBinStateEmissions(indsObsFirstState{:}))
%             % OK, 1st time we access that element, we calculate it
%             observation = data.(chrList{index}).binData(index2,:);
%             combinationsBinStateEmissions(indsObsAllStates{:}) =...
%                 getlogObsEmission(observation, logEmissions, logOneMinusEmissions);
%         end
%         logEmissionsBin.(chrList{index})(:,index2) = ...
%             squeeze(combinationsBinStateEmissions(indsObsAllStates{:}));
%     end
% end
% 
% modelOut.logEmissionsBin.(chrList{index}) = logEmissionsBin.(chrList{index});
% 

% 
% 
% function logObsEmission = getlogObsEmission(observation, logEmissions, logOneMinusEmissions)
% % Helper function to calculate the log of the probability of emmitting an
% % observation, given an emission matrix
% % Useful because there are a finite number of observations: 2^nMarks
% logObsEmission = sum(logEmissions.*repmat(observation == 1,size(logEmissions,1),1) +...
%     logOneMinusEmissions.*repmat(observation == 0,size(logEmissions,1),1),2);


function nTransitions = getNTransitions(listStates,state1,state2)
% Given a list and two states calculate the number of transitions

hit1States = listStates==state1;
if  state1==state2
    nTransitions = sum((2*hit1States(2:end)-hit1States(1:end-1))==1);
else
    hit2States = listStates==state2;
    nTransitions = sum((2*hit2States(2:end)-hit1States(1:end-1))==1);
end










