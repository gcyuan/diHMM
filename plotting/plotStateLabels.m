function plotStateLabels(modelIn, orientation, indices, stateLabels, labelType)
%   PLOTSTATELABELS helper function to display nice labels for states
%
%   Author: Eugenio Marco


if ~exist('orientation', 'var')
    orientation = 'horizontal';
end

if ~exist('stateLabels', 'var')
    nB = modelIn.param.nB;
    stateLabels = cell(nB,1);
    for index = 1:nB
        stateLabels{index} = ['N' num2str(index)];
    end
end

if ~exist('labelType', 'var')
    labelType = 'bin';
end

if ~exist('indices', 'var')
    indices = 1:modelIn.param.nB;
end


switch labelType
    case 'bin'
        allAnnotations = getBinAnnotations;
        if isfield(modelIn, 'assignedBinAnnotations')
            [statesAnnotations, dummy] = sort(modelIn.assignedBinAnnotations);
        else
            statesAnnotations = ones(1,modelIn.param.nB)*14; % 14 is WHITE
        end
    case 'domain'
        allAnnotations = getDomainAnnotations;
        if isfield(modelIn, 'assignedDomainAnnotations')
            [statesAnnotations, dummy] = sort(modelIn.assignedDomainAnnotations);
        else
            statesAnnotations = ones(1,modelIn.param.nD)*14; % 14 is WHITE
        end
end

annotationsColor = reshape([allAnnotations.color],3,length(allAnnotations))'/255;


%%
axes
switch orientation
    case 'vertical'
        imagesc(statesAnnotations(indices)',[1 size(annotationsColor,1)])
        set(gca,'yticklabel',stateLabels(indices),'xticklabel',[],'ytick',1:length(indices))
        set(gca,'pos',[0.1300    0.1100    0.0575    0.8150])
        colormap(annotationsColor)
        doLabels(stateLabels(indices))
        set(gca,'yticklabel',[],'xticklabel',[], 'ytick',[], 'xtick', [])
        freezeColors
    case'horizontal'
        imagesc(statesAnnotations(indices),[1 size(annotationsColor,1)])
        set(gca,'xticklabel',stateLabels(indices),'yticklabel',[],'xtick',1:length(indices))
        set(gca,'pos',[0.1300    0.1100    0.7750    0.0733])
        colormap(annotationsColor)
        doLabelsHorizontal(stateLabels(indices))
        set(gca,'yticklabel',[],'xticklabel',[], 'ytick',[], 'xtick', [])
        freezeColors
end

%%

end


function doLabels(stateLabelsR)

nD = length(stateLabelsR);
%erase current tick labels from figure
set(gca,'YTickLabel',[]);
%get tick label positions
b=get(gca,'XTick');
c=get(gca,'YTick');
%make new tick labels
text(repmat(b(1)+1*(b(2)-b(1)),nD,1),c,stateLabelsR,'HorizontalAlignment','center','fontsize',14);

end


function doLabelsHorizontal(stateLabelsR)

nD = length(stateLabelsR);
%erase current tick labels from figure
set(gca,'XTickLabel',[]);
%get tick label positions
c=get(gca,'XTick');
b=get(gca,'YTick');
%make new tick labels
text(c,repmat(b(1)+1*(b(2)-b(1)),1,nD),stateLabelsR','HorizontalAlignment','center','fontsize',14,'rotation',90);

end