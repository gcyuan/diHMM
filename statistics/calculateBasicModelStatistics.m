function calculateBasicModelStatistics(cellTypes, nB, nD, projectName, runNumber)
%   CALCULATEBASICMODELSTATISTICS Function that creates plots for state coverages and 
%   enrichments of domains into nucleosome-level states
%
%   Author: Eugenio Marco

baseName = getRunBaseName(projectName, cellTypes, runNumber, nB, nD);
load([baseName '.mat'])

% reorder the model
[~, newBinOrder] = sort(modelFinal.assignedBinAnnotations);
[~, newDomainOrder] = sort(modelFinal.assignedDomainAnnotations);
modelFinal = reorderModel(modelFinal, newBinOrder, newDomainOrder);

% get figures path and create directory if not present
figuresBaseName = getResultsFileBaseName(projectName, 'figures', [], runNumber, nB, nD);
[pathFigures, ~, ~] = fileparts(figuresBaseName);
checkBasename(figuresBaseName)

binLabels = cell(nB,1);
for index = 1:nB
    binLabels{index} = ['N' num2str(index)];
end


domainLabels = cell(nD,1);
for index = 1:nD
    domainLabels{index} = ['D' num2str(index)];
end

chrList = fieldnames(modelFinal.states.(cellTypes{1}));

binSize = modelFinal.param.binSize;
cmBlues = cbrewer('seq', 'Blues',100);
cmGreens = cbrewer('seq', 'Greens',100);
enrichmentScaleMax = 50;
enrichmentScaleMin = 1/enrichmentScaleMax;

nCT = length(cellTypes);

sizeBinLevelStates = cell(nCT,1);
sizeDomainLevelStates = cell(nCT,1);
domainLevelComposition = cell(nCT,1);
genomeSize= cell(nCT,1);
nDSize = modelFinal.param.nDSize;

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

%% Coverage Image Domains

allDomainSizes = [];
for index = 1:nCT
    allDomainSizes = [allDomainSizes sizeDomainLevelStates{index}/genomeSize{1}];
end

H1 = figure('position',[ 889    69   425   718], 'Color', 'w');
% axes('pos',[0.2500    0.0378    0.4949    0.8690]) % long
% HA = axes('pos',[0.3153    0.1100    0.4266    0.6630]);
HA = axes('pos',[0.3153    0.2200    0.3201    0.6630]);
nchalf = 4;
imagesc(log10(allDomainSizes), [-nchalf 0]);
cmBlues8 = cbrewer('seq', 'Blues',2*nchalf);
colormap(cmBlues8)
HT = text(1:nCT, nD*1.16*ones(1,3),...
    {cellTypes{2}, cellTypes{1}, cellTypes{3}},...
    'HorizontalAlignment','right','fontsize',18,'rotation',90,...
    'interpreter','none');
text(nCT/2.5,   .4, 'Coverage',...
    'HorizontalAlignment','left','fontsize',18,'rotation',00,...
    'interpreter','none');

freezeColors
HCB = colorbar;
cbfreeze(HCB)
ax = findobj(gcf,'Type','axes');


set(ax(1),'yticklabel',[], 'xticklabel', []) %Remove tick labels
set(HA,'ytick',[],'xtick',[])
% Get tick mark positions
yTicks = get(ax(1),'ytick');
xPos = 4.3936;
yPos = (0:(2*nchalf))/(2*nchalf)*nD+0.5128;
% Reset the ytick labels in desired font
niceLabels = {'10^{0}', ' ', '10^{-1}', ' ', '10^{-2}', ' ', '10^{-3}', ' ', '10^{-4}', ' ', '10^{-5}', ' ', '10^{-6}'};
for i = 1:length(yPos)
    %Create text box and set appropriate properties
    text(xPos,yPos(i),['$' niceLabels{i} '$'],...
        'HorizontalAlignment','Left','interpreter', 'latex', 'fontsize',18,...
        'FontName', 'Helvetica');
end
set(ax(1),'pos', [0.6624    0.22    0.0558    0.6630])

set(HA,'pos', [0.3153    0.2200    0.3201    0.6630])
plotStateLabels(modelFinal,'vertical', 1:nD, domainLabels, 'domain')
freezeColors
set(gca,'pos', [0.2153    0.22    0.1000    0.6630])

filenameeps = fullfile(pathFigures, [projectName '_' 'combined_coverages_domains_nB'...
    num2str(nB) '_nD' num2str(nD) '.eps']);
export_fig(filenameeps)


%% Coverage Image Bins

allBinSizes = [];
for index = 1:nCT
    allBinSizes = [allBinSizes sizeBinLevelStates{index}/genomeSize{1}];
end

H1 = figure('position',[ 889    69   425   718], 'Color', 'w');
% axes('pos',[0.2500    0.0378    0.4949    0.8690]) % long
% HA = axes('pos',[0.3153    0.1100    0.4266    0.6630]);
HA = axes('pos',[0.3153    0.2200    0.3201    0.6630]);
nchalf = 4;
imagesc(log10(allBinSizes), [-nchalf 0]);
cmBlues8 = cbrewer('seq', 'Blues',2*nchalf);
colormap(cmBlues8)
    
HT = text(1:nCT, nB*1.16*ones(1,3),...
    cellTypes,...
    'HorizontalAlignment','right','fontsize',18,'rotation',90,...
    'interpreter','none');
text(0.8383,   -0.6302, 'Coverage',...
    'HorizontalAlignment','left','fontsize',18,'rotation',00,...
    'interpreter','none');

freezeColors
HCB = colorbar;
cbfreeze(HCB)
ax = findobj(gcf,'Type','axes');


set(ax(1),'yticklabel',[], 'xticklabel', []) %Remove tick labels
set(HA,'ytick',[],'xtick',[])
% Get tick mark positions
yTicks = get(ax(1),'ytick');
xPos = 4.3936;
yPos = (0:(2*nchalf))/(2*nchalf)*nB+0.5128;
% Reset the ytick labels in desired font
niceLabels = {'10^{0}', ' ', '10^{-1}', ' ', '10^{-2}', ' ', '10^{-3}', ' ', '10^{-4}', ' ', '10^{-5}', ' ', '10^{-6}'};
for i = 1:length(yPos)
    %Create text box and set appropriate properties
    text(xPos,yPos(i),['$' niceLabels{i} '$'],...
        'HorizontalAlignment','Left','interpreter', 'latex', 'fontsize',18,...
        'FontName', 'Helvetica');
end
set(ax(1),'pos', [0.6624    0.22    0.0558    0.6630])

set(HA,'pos', [0.3153    0.2200    0.3201    0.6630])
plotStateLabels(modelFinal,'vertical', 1:nB, binLabels, 'bin')
freezeColors
set(gca,'pos', [0.2153    0.22    0.1000    0.6630])

filenameeps = fullfile(pathFigures, [projectName '_' 'combined_coverages_bins_nB'...
    num2str(nB) '_nD' num2str(nD) '.eps']);
export_fig(filenameeps)


%%  Combined Nucleosome-level composition

dataDLC = zeros(nD,nB);

for indexCT = 1:nCT
    dataDLC = dataDLC + domainLevelComposition{indexCT};
end

for indexD = 1:nD
    dataDLC(indexD,:) = dataDLC(indexD,:)/sum(dataDLC(indexD,:));
end

%% Domain enrichments


domainBinEnrichments = getDomainBinEnrichments(modelFinal);


%% Domain log10 fold enrichments masked by FDR

dataPValues = domainBinEnrichments.allCombinedDomainBinEnrichmentsPValues;


enrichmentFDRMaskedScaleMax = 50;
% enrichmentTimeslogPValueScaleMin = enrichmentTimeslogPValueScaleMax/100;


% dataCombinedPValues = -log10(dataPValues);
% dataCombinedPValues(isinf(dataCombinedPValues)) = 300; % Max p-value observed

% pValueThreshold = 200;
FDRThreshold = 0.01;

dataFDR =reshape(mafdr(dataPValues(:)),size(dataPValues));

dataEnrich = domainBinEnrichments.allCombinedDomainBinEnrichments;
dataEnrichmentFDRMasked = dataEnrich;
dataEnrichmentFDRMasked(dataFDR > FDRThreshold) = 0;

H1 = figure('position',[ 39          80        1415         718],'Color', 'w');
xCoord = 0.0571;
HA = axes('pos',[xCoord 0.18 0.82 0.73]);
imagesc(log10(dataEnrichmentFDRMasked+eps),[log10(enrichmentFDRMaskedScaleMax/100) log10(enrichmentFDRMaskedScaleMax)])

colormap(cmGreens)
    
title('Combined Nucleosome-level Fold Enrichments for each Domain','fontsize',18,'interpreter','none')
set(gca,'ytick',1:nD,'yticklabel',domainLabels)

set(gca,'fontsize',18)
set(gca,'xtick',1:nB,'xticklabel',binLabels)
HCB = colorbar;
cbfreeze(HCB)
freezeColors
ax = findobj(gcf,'Type','axes');


set(ax(1), 'xticklabel', []) %Remove tick labels
set(ax(1),'ytick', [log10(enrichmentFDRMaskedScaleMax/100) log10(enrichmentFDRMaskedScaleMax/10) log10(enrichmentFDRMaskedScaleMax)],...
    'yticklabel', [enrichmentFDRMaskedScaleMax/100 enrichmentFDRMaskedScaleMax/10 enrichmentFDRMaskedScaleMax],'fontsize',18) %Remove tick labels
set(HA,'ytick',[],'xtick',[])
% Get tick mark positions
yTicks = get(ax(1),'ytick');
xPos = 32.05;
yPos = nD:-(nD/4):-1+0.5128;

set(ax(1),'pos', [0.8456    0.1776    0.0168    0.7312])

axPos = get(ax(end),'pos');
newWidth = axPos(3);
set(gca,'ytick',[],'Xtick',[])
freezeColors
plotStateLabels(modelFinal,'horizontal', 1:nB, binLabels, 'bin')
set(gca,'pos',[xCoord    0.1100    newWidth    0.0733])
freezeColors
axes(HA)
plotStateLabels(modelFinal,'vertical', 1:nD, domainLabels, 'domain')
set(gca,'pos',[0.0233    0.1800    0.0320    0.7300])


filenameeps = fullfile(pathFigures, [projectName '_' 'Combined_Domain_Nucleosome_log_enrichments_masked_FDR_nB' num2str(nB) '_nD' num2str(nD) '.eps']);
export_fig(filenameeps)




%%
% Nice emissions

%%
if strcmp(projectName,'modelGHKChr17Final');
    reordering = [7, 6, 5, 8, 2, 4, 9, 1, 3];
else
    reordering = 1:length(modelFinal.markNames);
end
% cmBlues = cbrewer('seq', 'Blues',10);
cmPurples = cbrewer('seq', 'Purples',10);

H1 = figure('position',[477    58   660   754],'Color', 'w');
% axes('pos',[0.2500    0.0378    0.4949    0.8690]) % long
axes('pos',[0.2500    0.1631    0.65    0.7437])
imagesc(modelFinal.emissions(:,reordering), [0 1]);
nMarks = length(modelFinal.markNames);
% text(3.5, -1,'Emissions','fontsize',18)
set(gca,'ytick',1:nB,'Xtick',1:nMarks)
c=get(gca,'XTick');
b=get(gca,'YTick');
set(gca, 'xtick', [],'xticklabel', [])
%     set(gca,'xtick',1:nAnnotations)
%make new tick labels
HT =     text(c,repmat(b(end)*1.04,1,nMarks),modelFinal.markNames(reordering)',...
    'HorizontalAlignment','right','fontsize',14,'rotation',90,...
    'interpreter','none');
set(gca,'yticklabel',[],'xticklabel',[])
set(gca,'fontsize',18)
colormap(cmPurples)
set(gca,'ytick',[],'Xtick',[])
freezeColors
HCB = colorbar;
cbfreeze(HCB)
title('Nucleosome-level emissions')
plotStateLabels(modelFinal,'vertical', 1:nB, binLabels, 'bin')
set(gca,'pos',[0.1757    0.1631    0.0739    0.7437],'fontsize',18)
freezeColors

filenameeps = fullfile(pathFigures,[projectName '_' 'emissions_nB' num2str(nB) '_nD' num2str(nD) '.eps']);
export_fig(filenameeps)

%%

%
domainEmissions = cell(nCT,1);
dataDLC = cell(nCT,1);
dataDLCCopy = cell(nCT,1);

cmBlues = cbrewer('seq', 'Blues',10);
for indexCT = 1:nCT
    dataDLC{indexCT} = zeros(nD,nB);
    for indexD = 1:nD
        dataDLC{indexCT}(indexD,:) = domainLevelComposition{indexCT}(indexD,:)/sum(domainLevelComposition{indexCT}(indexD,:));
    end
    dataDLCCopy(indexCT) = dataDLC(indexCT);
end

for indexCT = 1:nCT
    H1 = figure('position',[477    58   760   754],'Color', 'w');
    axes('pos',[0.2500    0.1631    0.65    0.7437])
    domainEmissions{indexCT} = dataDLCCopy{indexCT}*modelFinal.emissions;
    imagesc(domainEmissions{indexCT}, [0 1]);
    colormap(cmBlues)
    set(gca,'ytick',1:nD,'Xtick',1:nMarks)
    c=get(gca,'XTick');
    b=get(gca,'YTick');
    set(gca, 'xtick', 1:nMarks,'xticklabel', [])
    %make new tick labels
    HT =     text(c,repmat(b(end)*1.14,1,nMarks),modelFinal.markNames',...
        'HorizontalAlignment','right','fontsize',14,'rotation',90,...
        'interpreter','none');
    set(gca,'yticklabel',domainLabels,'xticklabel',[])
    %     doDomainLabelsLeft(domainLabels)
    set(gca,'ytick',[],'Xtick',[])
    set(gca,'fontsize',18)
    freezeColors
    HCB = colorbar;
    cbfreeze(HCB)
    title(['Mark Composition for each Domain in ' cellTypes{indexCT}],'interpreter','none');
    
    plotStateLabels(modelFinal,'vertical', 1:nD, domainLabels, 'domain')
    set(gca,'pos',[0.1829    0.1631    0.0606    0.7437])
    
    
    filenameeps = fullfile(pathFigures, [projectName '_' 'mark_composition_' cellTypes{indexCT} '_nB' num2str(nB) '_nD' num2str(nD) '.eps']);
    export_fig(filenameeps)
    
end
%%

computeAllSizes(modelFinal,projectName, pathFigures)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function computeAllSizes(model,projectName, pathFigures)

if model.param.nD > 1
    computeSizes(model, 'domainLevel',projectName, pathFigures)
end

computeSizes(model, 'binLevel',projectName, pathFigures)

end
%%%%%%%
function computeSizes(model, level, projectName, pathFigures)


nDSize = model.nDSize;
binSize = model.binSize;
cellTypes = fields(model.states);
nCT = length(cellTypes);

chrList = fieldnames(model.states.(cellTypes{1}));


for indexCT = 1:nCT
    
    
    for index = 1:length(chrList)
        switch level
            case 'binLevel'
                states = model.states.(cellTypes{indexCT}).(chrList{index}).binStates;
                upperBinLocations = binSize:binSize:size(states,1)*binSize';
            case 'domainLevel'
                states = model.states.(cellTypes{indexCT}).(chrList{index}).domainStates;
                upperBinLocations = binSize*nDSize:nDSize*binSize:nDSize*size(states,1)*binSize';
        end
        
        transitions = (double(states(1:end-1))-double(states(2:end))~=0);
        indexTransitions = find(transitions);
        
        upperStateLocations = [upperBinLocations(indexTransitions) upperBinLocations(end)]';
        bottomStateLocations = [0 upperStateLocations(1:end-1)']';
        stateNames = [states(transitions); states(indexTransitions(end)+1)];
        stateSizes = upperStateLocations - bottomStateLocations;
        
        if index == 1
            allStateNames = stateNames;
            allStateSizes = stateSizes;
        else
            allStateNames = [allStateNames; stateNames];
            allStateSizes = [allStateSizes; stateSizes];
        end
        
    end
    %
    
    
    
    if indexCT == 1
        allCombinedSizesNames = [allStateSizes, allStateNames];
    else
        allCombinedSizesNames = [allCombinedSizesNames; [allStateSizes, allStateNames]];
    end
    

    nStates = size(unique(allStateNames),1);
    stateMedianSizes = zeros(nStates,2);
    for index = 1:nStates
        stateMedianSizes(index,:) = [index,median(allStateSizes(allStateNames==index))];
        
    end
    disp(['stateMedianSizes for ' level])
    disp(' ')
    disp(stateMedianSizes)
    
    stateMeanSizes = zeros(nStates,2);
    for index = 1:nStates
        stateMeanSizes(index,:) = [index,mean(allStateSizes(allStateNames==index))];
        
    end
    disp(' ')
    
    disp(['stateMeanSizes for ' level])
    disp(' ')
    disp(stateMeanSizes)
    disp(' ')
end


%%
H1 = figure('position',[ 39          80        1415         718],'Color', 'w');

boxplot(log10(allCombinedSizesNames(:,1)),allCombinedSizesNames(:,2))
set(gca,'xtick',[], 'fontsize', 18)
ylabel('log10(Length) (bases)')
set(gca, 'ytick',2:7)

switch level
    case 'binLevel'
        ylim([2 6])
        hold on
        nB = model.param.nB;
%         yValue = log10(model.param.binSize);
%         plot([0 nB+1], [yValue yValue],'--')
        binLabels = cell(nB,1);
        for index = 1:nB
            binLabels{index} = ['N' num2str(index)];
        end
        plotStateLabels(model,'horizontal', 1:nB, binLabels, 'bin')
    case 'domainLevel'
        ylim([3.5 6.6])
        nD = model.param.nD;
        hold on
%         yValue = log10(model.param.nDSize* model.param.binSize);
%         plot([0 nD+1], [yValue yValue],'k--')
        domainLabels = cell(nD,1);
        for index = 1:nD
            domainLabels{index} = ['D' num2str(index)];
        end
        plotStateLabels(model,'horizontal', 1:model.param.nD, domainLabels, 'domain')
end

set(gca,'pos', [0.1300    0.0571    0.7750    0.0733])
freezeColors


filenameeps = fullfile(pathFigures, [projectName '_' 'sizes_' level '_' num2str(model.param.nB) '_nD' ...
    num2str(model.param.nD) '.eps']);
export_fig(filenameeps)


%%


end

