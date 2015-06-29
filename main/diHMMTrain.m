function results = diHMMTrain(data, initialModel, param, chrList)
%   DIHMMTRAIN function to train the Hierarchical Hidden Markov Model, diHMM
%   Calculates a results structure containing the scaled forward and backward
%   variables, scaling factors, posterior probabilities and loglikelihood for diHMM
%
%   Author: Eugenio Marco

tol = param.tol;
maxIter = param.maxIter;

nB = param.nB;
nD = param.nD;
tB = initialModel.transitionB;
tD = initialModel.transitionD;

cellTypes = fields(data)';
nCellTypes = length(fields(data));

% Preallocate
initialModel.logLik = 0; %Dummy value
submodel(maxIter+1) = extractParameters(initialModel);
submodel(maxIter+1) = [];

runBaseName = getRunBaseName(param.projectName,cellTypes, param.runNumber, nB, nD);
checkBasename(runBaseName)


% em = initialModel.emissionsBin2dec;

% tDExpanded = kron(tD, ones(nB));
% tDIntraDomain = kron(diag(ones(nD,1)), ones(nB));

% tBExpanded = repmat(reshape(tB,nB,nB*nD),nD,1);

estimatedTB = zeros(size(tB));
estimatedTD = zeros(size(tD));
estimatedEmissions = zeros(size(initialModel.emissions));

% INDEX ORDER: BinState, Domain, BinPosition
% forwardScaled = zeros(nB,nD,nBinPosition);
% backwardScaled = zeros(nB,nD,nBinPosition);
% scale = zeros(1,nBinPosition);


converged = false;
logLik = 1;
logLiks = zeros(1,maxIter);

updatedModel = initialModel;

nChr = length(chrList);

% logLikCTChr = zeros(nCellTypes,nChr,1);

for iteration = 1:maxIter
    %     numberTransitionsTB = zeros(size(estimatedTB));
    %     numberTransitionsTD = zeros(size(estimatedTD));
%     if nD==1
%         % We need an extra singleton
%         numberTransitionsTBCTChr = zeros([size(estimatedTB),1,nCellTypes,nChr]);
%     else
%         numberTransitionsTBCTChr = zeros([size(estimatedTB),nCellTypes,nChr]);
%     end
%     numberTransitionsTDCTChr = zeros([size(estimatedTD),nCellTypes,nChr]);
%     numberEmissionsCTChr = zeros([size(estimatedEmissions),nCellTypes,nChr]);
%     expectedBinStatesCTChr = zeros(nB,nCellTypes,nChr);
%     estimatedInitialProbabilitiesCTChr = zeros(nB,nD,nCellTypes,nChr);
    
    % We linearize the indices for CT and Chr, to run it fully parallel
    if nD==1
        % We need an extra singleton
        numberTransitionsTBCTChrIntermediate = zeros([size(estimatedTB),1,nCellTypes*nChr]);
    else
        numberTransitionsTBCTChrIntermediate = zeros([size(estimatedTB),nCellTypes*nChr]);
    end
    numberTransitionsTDCTChrIntermediate = zeros([size(estimatedTD),nCellTypes*nChr]);
    numberEmissionsCTChrIntermediate = zeros([size(estimatedEmissions),nCellTypes*nChr]);
    expectedBinStatesCTChrIntermediate = zeros(nB,nCellTypes*nChr);
    estimatedInitialProbabilitiesCTChrIntermediate = zeros(nB,nD,nCellTypes*nChr);
    logLikCTChrIntermediate = zeros(nCellTypes*nChr,1);
    
    oldLL = logLik;
    
    updatedModel.tDExpanded = kron(updatedModel.transitionD, ones(nB));
    
    updatedModel.tBExpanded = repmat(reshape(updatedModel.transitionB,nB,nB*nD),nD,1);
    
    parfor indexCTChr = 1:nCellTypes*nChr
        % We run in parallel all chromosomes
        [numberTransitionsTBCTChrIntermediate(:,:,:,indexCTChr), numberTransitionsTDCTChrIntermediate(:,:,indexCTChr),...
            numberEmissionsCTChrIntermediate(:,:,indexCTChr), expectedBinStatesCTChrIntermediate(:,indexCTChr),...
            estimatedInitialProbabilitiesCTChrIntermediate(:,:,indexCTChr), logLikCTChrIntermediate(indexCTChr)] = ...
            runCTChr(data, updatedModel, param, indexCTChr);
    end
    
    %   Now we unfold the matrices into their indexCT and indexChr parts
    numberTransitionsTBCTChr = reshape(numberTransitionsTBCTChrIntermediate, ...
        [size(estimatedTB),nCellTypes,nChr]);
    numberTransitionsTDCTChr = reshape(numberTransitionsTDCTChrIntermediate, ...
        [size(estimatedTD),nCellTypes,nChr]);
    numberEmissionsCTChr = reshape(numberEmissionsCTChrIntermediate, ...
        [size(estimatedEmissions),nCellTypes,nChr]);
    expectedBinStatesCTChr = reshape(expectedBinStatesCTChrIntermediate, nB,nCellTypes,nChr);
    estimatedInitialProbabilitiesCTChr = reshape(estimatedInitialProbabilitiesCTChrIntermediate, ...
        nB,nD,nCellTypes,nChr);
    logLikCTChr = reshape(logLikCTChrIntermediate, nCellTypes,nChr,1);
    
    
    
    numberTransitionsTB = sum(sum(numberTransitionsTBCTChr,5),4);
    numberTransitionsTD = sum(sum(numberTransitionsTDCTChr,4),3);
    numberEmissions = sum(sum(numberEmissionsCTChr,4),3);
    expectedBinStates = sum(sum(expectedBinStatesCTChr,3),2);
    estimatedInitialProbabilities = sum(sum(estimatedInitialProbabilitiesCTChr,4),3)/nCellTypes/nChr;
    logLik = sum(logLikCTChr(:));
    
    % Reestimate tB
    totalNumberTransitionsTB = sum(numberTransitionsTB,2);
    estimatedTransitionB = numberTransitionsTB./totalNumberTransitionsTB(:,ones(nB,1),:);
    for indexD = 1:param.nD
        if any(totalNumberTransitionsTB(:,1,indexD) == 0)
            matrixFormTB = estimatedTransitionB(:,:,indexD);
            noTransitionRows = find(totalNumberTransitionsTB(:,1,indexD) == 0);
            matrixFormTB(noTransitionRows,:) = 0;
            matrixFormTB(sub2ind(size(matrixFormTB),noTransitionRows,noTransitionRows)) = 1;
            estimatedTransitionB(:,:,indexD) = matrixFormTB;
        end
    end
    % clean up any remaining Nans
    estimatedTransitionB(isnan(estimatedTransitionB)) = 0;
    
    % Reestimate tD
    totalNumberTransitionsTD = sum(numberTransitionsTD,2);
    estimatedTransitionD = numberTransitionsTD./repmat(totalNumberTransitionsTD,1,nD);
    if any(totalNumberTransitionsTD == 0)
        noTransitionRows = find(totalNumberTransitionsTD == 0);
        estimatedTransitionD(noTransitionRows,:) = 0;
        estimatedTransitionD(sub2ind(size(estimatedTransitionD),noTransitionRows,noTransitionRows)) = 1;
    end
    % clean up any remaining Nans
    estimatedTransitionD(isnan(estimatedTransitionD)) = 0;
    
    
    % Reestimate emissions
    estimatedEmissions = numberEmissions./repmat(expectedBinStates,1,param.nMarks);
    
    
    updatedModel.emissions = estimatedEmissions;
    updatedModel.transitionB = estimatedTransitionB;
    updatedModel.transitionD = estimatedTransitionD;
    updatedModel.initialProbabilities = estimatedInitialProbabilities;
    updatedModel.logLik = logLik;
    
    % Update logTransitionB and so on, but check later if needed
    
    updatedModel = postProcessModel(updatedModel, param);
    logLiks(iteration) = logLik;

    updatedModel.logLiks = logLiks;
    updatedModel.iteration = iteration;
    
    save([runBaseName 'Intermediate.mat'], 'updatedModel')
    
    submodel(iteration) = extractParameters(updatedModel);
    save([runBaseName 'Submodel.mat'], 'submodel')
    
        
    if (abs(logLik-oldLL)/(1+abs(oldLL))) < tol
        % We keep this in case we want to implement it
        %         if norm(guessTR - oldGuessTR,inf)/numStates < trtol
        %             if norm(guessE - oldGuessE,inf)/numEmissions < etol
        %                 if verbose
        %                     fprintf('%s\n',getString(message('stats:hmmtrain:ConvergedAfterIterations',iteration)))
        %                 end
        converged = true;
        break
        %             end
        %         end
    end
    
    if iteration > 1
        disp(['Iteration = ' num2str(iteration) '  LogLikelihood = ' num2str(logLiks(iteration)) '  Change = ' ...
            num2str(logLiks(iteration)-logLiks(iteration-1))])
    else
        disp(['LogLikelihood = ' num2str(logLiks(iteration))])
    end
end


if ~converged
    warning('No Convergence');
end

results = updatedModel;
end


function [numberTransitionsTBCTChrIntermediate, numberTransitionsTDCTChrIntermediate,...
    numberEmissionsCTChrIntermediate, expectedBinStatesCTChrIntermediate,...
    estimatedInitialProbabilitiesCTChrIntermediate, logLikCTChrIntermediate]  = ...
    runCTChr(data, updatedModel, param, indexCTChr)

nCT = length(param.cellTypes);
chrList = fields(data.(param.cellTypes{1}));
% nChr = length(chrList);

% We decompose the index indexCTChr into (indeCT, indexChr) pair
% indexCT = mod(indexCTChr+1, nCT) +1;
% indexChr = mod(indexCTChr-nCT*(+1, nChr) +1;
indexCT = rem(indexCTChr+1, nCT) +1;
indexChr = floor((indexCTChr-1)/nCT) +1;

dataCT = data.(param.cellTypes{indexCT});

[numberTransitionsTBCTChrIntermediate, numberTransitionsTDCTChrIntermediate,...
    numberEmissionsCTChrIntermediate, expectedBinStatesCTChrIntermediate,...
    estimatedInitialProbabilitiesCTChrIntermediate, logLikCTChrIntermediate] = ...
    runChr(dataCT, updatedModel, param, chrList{indexChr});

end


function [numberTransitionsTB, numberTransitionsTD,...
    numberEmissions, expectedBinStates,...
    estimatedInitialProbabilities, logLik]  = runChr(data, model, param, chr)

nB = param.nB;
nD = param.nD;

numberTransitionsTB = zeros(size(model.transitionB));
numberTransitionsTD = zeros(size(model.transitionD));
numberEmissions = zeros(size(model.emissions));

nBinPosition = size(data.(chr).binData,1);
bin2decData = data.(chr).bin2decData;
%     transition =

decoding = diHMMDecode(data, model, param, chr);
%         logf = log(decoding.forwardScaled);
%         logb = log(decoding.backwardScaled);
logLik = decoding.logLikelihood;

% We now need to calculate for each position the probability of
% transitioning from Bin = j, Dom = mu, to Bin = k, Dom = nu
% To vectorize we work with a vector, instead of (Bin,Dom) coordinates

for indexBinPosition = 1:nBinPosition-1
    
    if (mod(indexBinPosition, param.nDSize)~=0) || nD == 1
        % Intradomain case
        for indexD = 1:nD
            numberTransitionsTB(:,:,indexD) = numberTransitionsTB(:,:,indexD) +...
                decoding.forwardScaled(:,indexD*ones(1,nB),indexBinPosition) .* ...
                model.transitionB(:,:,indexD) .* ...
                model.emissionsBin2dec(:,(bin2decData(indexBinPosition+1)+1)*ones(1,nB))' .* ...
                decoding.backwardScaled(:,indexD*ones(1,nB),indexBinPosition+1)'/...
                decoding.scale(indexBinPosition+1);
        end
    else % Interdomain case
        boundaryTransitionMatrix = ...
            reshape(decoding.forwardScaled(:,:,indexBinPosition,ones(1,nB*nD)), nB*nD,nB*nD) .* ...
            model.tDExpanded .* model.tBExpanded .* ...
            reshape(model.emissionsBin2dec(:,bin2decData(indexBinPosition+1)+1,ones(1,nD*nD*nB)),nB*nD,nB*nD)' .* ...
            reshape(decoding.backwardScaled(:,:,indexBinPosition+1,ones(1,nB*nD)), nB*nD,nB*nD)'/...
            decoding.scale(indexBinPosition+1);
        reshapedBoundaryTransitionMatrix = reshape(boundaryTransitionMatrix, nB, nD, nB, nD);
        
        numberTransitionsTD = numberTransitionsTD +...
            squeeze(sum(sum(reshapedBoundaryTransitionMatrix,1),3));
        
        numberTransitionsTB = numberTransitionsTB +...
            squeeze(sum(reshapedBoundaryTransitionMatrix,2));
    end
end

% Reestimate emissions

posteriors = decoding.forwardScaled .* decoding.backwardScaled;
posteriorsBin = squeeze(sum(posteriors,2)); % Sum over tD

expectedBinStates = sum(posteriorsBin, 2); % Sum over all bins
for indexMark = 1:param.nMarks
    numberEmissions(:, indexMark) = sum(posteriorsBin(:,data.(chr).binData(:,indexMark)),2);
end

estimatedInitialProbabilities = posteriors(:,:,1);

end


function submodel = extractParameters(model)
submodel.emissions = model.emissions;
submodel.transitionB = model.transitionB;
submodel.transitionD = model.transitionD;
submodel.initialProbabilities = model.initialProbabilities;
submodel.logLik = model.logLik;
end


