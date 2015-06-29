function results = domHMMDecode(data, modelIn, param, chromosome)
% Calculates a results structure containing the scaled forward and backward
% variables, scaling factors, posterior probabilities and loglikelihood for domHMM


nBinPosition = size(data.(chromosome).binData,1);
nD = param.nD;
nB = param.nB;
bin2decData = data.(chromosome).bin2decData;
tB = modelIn.transitionB;
tD = modelIn.transitionD;
em = modelIn.emissionsBin2dec;

% INDEX ORDER: BinState, Domain, BinPosition
forwardScaled = zeros(nB,nD,nBinPosition);
backwardScaled = zeros(nB,nD,nBinPosition);
scale = zeros(1,nBinPosition);

% Initialization 
intermediateForwardScaled(:,:) = modelIn.initialProbabilities.*...
    reshape(em(:,bin2decData(1)+1,ones(1,nD)),nB,nD);
% We use reshape, instead of repmat. It's faster.
%     repmat(modelIn.emissionsBin2dec(:,bin2decData(1)+1),1,nD); % +1 since 0 is in element 1, 1 in element 2
scale(1) = sum(intermediateForwardScaled(:));
forwardScaled(:,:,1) = intermediateForwardScaled/scale(1);


% looping through binPositions
for indexBinPosition = 2:nBinPosition
    % We test whether we are at the domain boundary
    if (mod(indexBinPosition, param.nDSize)~=1) || nD == 1
        % Case within domain,
        % We loop through nD domains, for each we calculate a vector of
        % dimension nB
        for indexD=1:nD
            intermediateForwardScaled(:,indexD) = ...
                (forwardScaled(:, indexD, indexBinPosition-1)'*tB(:,:,indexD))';
        end
    else
        % Case of possible domain boundary
        for indexD=1:nD
            intermediateForwardScaled(:,indexD) = ...
                (sum(forwardScaled(:, :, indexBinPosition-1).*reshape(tD(:,indexD,ones(1,nB)),nD,nB)',2)'*...
                tB(:,:,indexD))';
        end
    end
    intermediateForwardScaled = intermediateForwardScaled.*...
        reshape(em(:,bin2decData(indexBinPosition)+1,ones(1,nD)),nB,nD);
    scale(indexBinPosition) = sum(intermediateForwardScaled(:));
    forwardScaled(:,:,indexBinPosition) = intermediateForwardScaled/scale(indexBinPosition);
end

% Begin backtracking
% Note that we shift one position the division by the scaling factor!!

% Initialization 
backwardScaled(:,:,nBinPosition) = 1;

for indexBinPosition = nBinPosition-1:-1:1
    % We test whether we are at the domain boundary
    if (mod(indexBinPosition, param.nDSize)~=0) || nD == 1
        % Case within domain,
        % We loop through nD domains, for each we calculate a vector of
        % dimension nB
        for indexD=1:nD
            intermediateBackwardScaled(:,indexD) = tB(:,:,indexD)*...
                (backwardScaled(:, indexD, indexBinPosition+1).*em(:,bin2decData(indexBinPosition+1)+1));
        end
    else
        % Case of possible domain boundary
        % Because of tB we need to separate now indexB instead of indexD
        for indexB=1:nB
            intermediateBackwardScaled(indexB,:) = ...
                (tD*...
                sum(backwardScaled(:, :, indexBinPosition+1).*...
                reshape(tB(indexB,:,:),nB,nD).*...
                reshape(em(:,bin2decData(indexBinPosition+1)+1,ones(1,nD)),nB,nD),1)')';
        end
        
        %         for indexD=1:nD
        %             for indexB=1:nB
        %                 intermediateBackwardScaled(indexB,indexD) = ...
        %                     sum(sum(backwardScaled(:, :, indexBinPosition+1).*...
        %                     reshape(em(:,bin2decData(indexBinPosition)+1,ones(1,nD)),nB,nD).*...
        %                     reshape(tB(indexB,:,:),nB,nD).*...
        %                     tD(indexD*ones(nB,1),:)));
        %             end
        %         end
    end
    backwardScaled(:, :, indexBinPosition) = ...
        intermediateBackwardScaled/scale(indexBinPosition+1); % SHIFT!!!
end

results.forwardScaled = forwardScaled;
results.backwardScaled = backwardScaled;
results.scale = scale;
results.logLikelihood = sum(log(scale));
    
 