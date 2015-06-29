function modelOut = postProcessModel(modelIn, param)
%   POSTPROCESSMODEL
%   Postprocession of the model matrices
%
%   Author: Eugenio Marco

nB = param.nB;
nM = length(modelIn.markNames);
% chrList = fieldnames(data);

modelIn.logEmissions = log([modelIn.emissions 1-modelIn.emissions]);
modelIn.logTransitionD = log(modelIn.transitionD);
modelIn.logTransitionB = log(modelIn.transitionB);

modelOut = modelIn;


% We store in dec format all emission probabilities
logEmissionsBin2dec = zeros(param.nB, 2^nM);
for index = 1:2^nM
    
    %OLD
%     % First column is equal to dec 0
%     binvec = fliplr(dec2bin(index-1,nM));
%     % From dec2binvec
%     % Convert the binary string, '1011', to a binvec, [1 1 0 1].
%     binvec = logical(str2num([fliplr(binvec);blanks(length(binvec))]')');

    %NEW
    % First column is equal to dec 0
    binvec = dec2bin(index-1,nM);
    % From dec2binvec
    % Convert the binary string, '1011', to a binvec, [1 1 0 1].
    binvec = logical(str2num([binvec;blanks(length(binvec))]')');
        

    intermediateLogEmissions = modelIn.logEmissions.*[...
        binvec(ones(nB,1),:), ~binvec(ones(nB,1),:)];
    % We eliminate the NaNs resulting from -Inf * 0
    intermediateLogEmissions(isnan(intermediateLogEmissions))=0;
    logEmissionsBin2dec(:,index) = sum(intermediateLogEmissions,2);
end
   
emissionsBin2dec = exp(logEmissionsBin2dec);

modelOut.logEmissionsBin2dec = logEmissionsBin2dec;
modelOut.emissionsBin2dec = emissionsBin2dec;
