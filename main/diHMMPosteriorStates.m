function modelOut = diHMMPosteriorStates(data, modelIn, param)
%   DIHMMPOSTERIORSTATES
%   Posterior State decoding for diHMM
%
%   Author: Eugenio Marco

modelOut = modelIn;
cellTypes = fields(data);
chrList = fieldnames(data.(cellTypes{1}));

for indexCT =1:length(cellTypes)
    dataCT = data.(cellTypes{indexCT});

    for index =1:length(chrList)
        % Each chromosome is independent. We find for each the best path
        % Following Rabiner we use phi for log(delta) and
        % define chi = argmax(phi +T)
        nBinPosition = size(dataCT.(chrList{index}).binData,1);
        nD = param.nD;
        nB = param.nB;
        
        decoding = diHMMDecode(dataCT, modelOut, param, chrList{index});
        
        posterior = decoding.forwardScaled .* decoding.backwardScaled;
        
        domainMarginal = squeeze(sum(posterior,1));
        posteriorsReshaped = reshape(posterior, nD*nB, nBinPosition);
        
        [maxPosterior, states] = max(posteriorsReshaped);
        
        binStates = rem(states-1,param.nB)+1;
        domainStatesSmallBins = floor((states-1)/param.nB)+1;
        
        modelOut.allLogLiks.(cellTypes{indexCT}).(chrList{index}) = decoding.logLikelihood;
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).binStates = binStates';
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainStatesSmallBins = domainStatesSmallBins';
        modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainStates = ...
            modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainStatesSmallBins(param.nDSize:param.nDSize:nBinPosition);
        if param.calculatePosteriorAndMarginal
            modelOut.states.(cellTypes{indexCT}).(chrList{index}).maxPosterior = maxPosterior;
            modelOut.states.(cellTypes{indexCT}).(chrList{index}).domainMarginal = domainMarginal;
        end
    end
end


