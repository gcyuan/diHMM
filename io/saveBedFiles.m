function saveBedFiles(baseName, model)
% SAVEBEDFILES function to generate bed files from the state calls
%
% Author: Eugenio Marco

[pathstr, ~, ~] = fileparts(baseName);

sepLocs = strfind(pathstr, filesep);
cellType = pathstr((sepLocs(end)+1):end);

filenameBinLevelStatesColor = [baseName  '_binLevelStatesColor.bed'];

filenameDomainLevelStatesColor = [baseName  '_domainLevelStatesColor.bed'];

filenameMaxPosterior = [baseName  '_posterior.bedgraph'];

saveFile(filenameBinLevelStatesColor, model, cellType, 'binLevel', 'color')

if model.param.nD > 1
    saveFile(filenameDomainLevelStatesColor, model, cellType, 'domainLevel', 'color')
end

if model.param.saveBedgraphPosterior
    saveBedgraphFile(filenameMaxPosterior, model, cellType, 'maxPosterior')
end

if model.param.saveBedgraphPosterior
    for index = 1:model.param.nD
        filenameDomainMarginal = [baseName  '_Domain_' num2str(index) '_marginal.bedgraph'];
        saveBedgraphFile(filenameDomainMarginal, model, cellType, 'domainMarginal', index)
    end
end


end

function saveFile(filename, model, cellType, level, varargin)

nDSize = model.param.nDSize;
binSize = model.param.binSize;

chrList = fieldnames(model.states.(cellType));

[~, name, ~] = fileparts(filename);


fs = fopen(filename, 'W+');

if nargin >4
    optionColor = varargin{1};
else
    optionColor = '';
end

trimLength = 100; % Length of domain trimming for visualization (in bp)

if strcmp(optionColor, 'color')
    fprintf(fs, ['track name=\"' name '\" description=\"' name...
        '\" visibility=1 itemRgb=\"On\"\n']);
    switch level
        case 'binLevel'
            binAnnotation = getBinAnnotations;
            annotationsColor = reshape([binAnnotation.color],3,length(binAnnotation))';
            if isfield(model,'assignedBinAnnotations') && ...
                    ~all(model.assignedBinAnnotations == model.assignedBinAnnotations(1))
                [statesAnnotations, ~] = sort(model.assignedBinAnnotations);
            else
                nColors = size(model.transitionB,1);
                colormap(hsv(nColors))
                cm = colormap;
                annotationsColor = uint8(255*cm);
                statesAnnotations = 1:model.param.nB;
            end
        case 'domainLevel'
            domainAnnotation = getDomainAnnotations;
            annotationsColor = reshape([domainAnnotation.color],3,length(domainAnnotation))';
            if isfield(model,'assignedDomainAnnotations')
                [statesAnnotations, ~] = sort(model.assignedDomainAnnotations);
            else
                nColors = size(model.transitionD,1);
                colormap(hsv(nColors))
                cm = colormap;
                annotationsColor = uint8(255*cm);
                statesAnnotations = 1:model.param.nD;
            end
    end
end

for index = 1:length(chrList)
    switch level
        case 'binLevel'
            states = model.states.(cellType).(chrList{index}).binStates;
            upperBinLocations = binSize:binSize:size(states,1)*binSize';
            stateSign = 'N';
        case 'domainLevel'
            states = model.states.(cellType).(chrList{index}).domainStates;
            upperBinLocations = binSize*nDSize:nDSize*binSize:nDSize*size(states,1)*binSize';
            stateSign = 'D';
    end
    
    transitions = (double(states(1:end-1))-double(states(2:end))~=0);
    indexTransitions = find(transitions);
    
    
    switch level
        case 'binLevel'
            upperStateLocations = [upperBinLocations(indexTransitions) upperBinLocations(end)]';
            bottomStateLocations = [0 upperStateLocations(1:end-1)']';
            stateNames = [states(transitions); states(indexTransitions(end)+1)];
         case 'domainLevel' % We trim domains
            upperStateLocations = [upperBinLocations(indexTransitions) upperBinLocations(end)]'-trimLength;
            bottomStateLocations = [0 (upperStateLocations(1:end-1)'++trimLength)]'+trimLength;
            stateNames = [states(transitions); states(indexTransitions(end)+1)];
    end
    
    for index2 = 1:length(upperStateLocations)
        fprintf(fs, '%s\t', chrList{index});
        fprintf(fs, '%d\t', bottomStateLocations(index2));
        fprintf(fs, '%d\t', upperStateLocations(index2));
        if strcmp(optionColor, 'color')
            fprintf(fs, '%s\t0\t.\t', [stateSign num2str(stateNames(index2))]);
            fprintf(fs, '%d\t', bottomStateLocations(index2));
            fprintf(fs, '%d\t', upperStateLocations(index2));
            fprintf(fs, '%d,%d,%d\n', annotationsColor(statesAnnotations(stateNames(index2)),:));
        else
            fprintf(fs, '%s\n', [stateSign num2str(stateNames(index2))]);
        end
    end
end

fclose all;

end


function saveIndividualDomainFile(baseName, model, cellType, splitDomains)

if ~exist('splitDomains','var')
    splitDomains = true;
end
nDSize = model.nDSize;
binSize = model.binSize;
nD = model.param.nD;

chrList = fieldnames(model.states.(cellType));

if nD > 1
    for indexND = 1:nD
        filename = [baseName '_individualDomain_' num2str(indexND) '_States.bed'];
        
        fs = fopen(filename, 'W+');
        
        for index = 1:length(chrList)
            
            states = model.states.(cellType).(chrList{index}).domainStates;
            upperBinLocations = binSize*nDSize:nDSize*binSize:nDSize*size(states,1)*binSize';
            stateSign = 'D';
            
            transitions = (double(states(1:end-1))-double(states(2:end))~=0);
            if splitDomains
                transitions(:) = true;
            end
            indexTransitions = find(transitions);
            
            upperStateLocations = [upperBinLocations(indexTransitions) upperBinLocations(end)]';
            bottomStateLocations = [0 upperStateLocations(1:end-1)']';
            stateNames = [states(transitions); states(indexTransitions(end)+1)];
            
            selectedND = stateNames==indexND;
            bottomStateLocations = bottomStateLocations(selectedND);
            upperStateLocations = upperStateLocations(selectedND);
            
            for index2 = 1:length(upperStateLocations)
                fprintf(fs, '%s\t', chrList{index});
                fprintf(fs, '%d\t', bottomStateLocations(index2));
                fprintf(fs, '%d\t', upperStateLocations(index2));
                fprintf(fs, '%s\n', [stateSign num2str(indexND)]);
            end
        end
    end
    fclose all;

end

end

function saveBedgraphFile(filename, model, cellType, type, domainLevel)

chrList = fieldnames(model.states.(cellType));
binSize = model.binSize;

switch type
    case 'maxPosterior'
        doBedgraph = isfield(model.states.(cellType).(chrList{1}), 'maxPosterior');
    case 'domainMarginal'
        doBedgraph = isfield(model.states.(cellType).(chrList{1}), 'domainMarginal');
end

if doBedgraph
    [~, name, ~] = fileparts(filename);
    
    
    fs = fopen(filename, 'W+');
    
    fprintf(fs, ['track type=bedGraph name=\"' name '\" description=\"' name...
        '\" visibility=full color=0,0,255 graphType=points\n']);
    
    
    for index = 1:length(chrList)
        switch type
            case 'maxPosterior'
                scores = model.states.(cellType).(chrList{index}).maxPosterior;
            case 'domainMarginal'
                scores = model.states.(cellType).(chrList{index}).domainMarginal(domainLevel,:);
        end
        upperBinLocations = (binSize:binSize:size(scores,2)*binSize)';
        bottomBinLocations = [0; upperBinLocations(1:end-1)];
        
        for index2 = 1:length(scores)
            fprintf(fs, '%s\t', chrList{index});
            fprintf(fs, '%d\t', bottomBinLocations(index2));
            fprintf(fs, '%d\t', upperBinLocations(index2));
            fprintf(fs, '%d\n', scores(index2));
        end
    end
    
    fclose all;
    
end

end

