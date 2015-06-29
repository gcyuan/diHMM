function plotDomainAnnotations
%   PLOTDOMAINANNOTATIONS function to display a legend with the available
%   domain annotations
%
%   Author: Eugenio Marco

domainAnnotation = getDomainAnnotations;

annotationsColor = reshape([domainAnnotation.color],3,length(domainAnnotation))'/255;

nAnnotations = size(annotationsColor,1);

annotationLabels = cell(nAnnotations,1);
for index = 1:nAnnotations
    annotationLabels{index} = domainAnnotation(index).name;
end

%%

figure
set(gcf,'pos',[520   378   441   420])
axes
imagesc((1:nAnnotations)',[1 nAnnotations])
set(gca,'yticklabel',1:nAnnotations,'xticklabel',[],'ytick',1:nAnnotations)
% set(gca,'pos',[0.1300    0.1100    0.0575    0.8150])
colormap(annotationsColor)
doLabels(annotationLabels)
set(gca,'yticklabel',[],'xticklabel',[], 'ytick',[], 'xtick', [])
freezeColors
   
%%

end


function doLabels(binLabelsR)

nD = length(binLabelsR);
%erase current tick labels from figure
set(gca,'YTickLabel',[]);
%get tick label positions
% b=get(gca,'XTick');
c=get(gca,'YTick');
%make new tick labels
text(repmat(.55,nD,1),c,binLabelsR,'HorizontalAlignment','left','fontsize',18);

end

