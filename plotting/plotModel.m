function plotModel(model)
%   PLOTMODEL plots model parameters given the model as input
%
%   Author: Eugenio Marco

% We reorder the model to displayed it sorted by annotation
[~, newBinOrder] = sort(model.assignedBinAnnotations);
[~, newDomainOrder] = sort(model.assignedDomainAnnotations);

model = reorderModel(model, newBinOrder, newDomainOrder);

param = model.param;

if isfield(model, 'markNames')
    markNames = model.markNames;
else
    markNames = param.markNames; %Hack for the moment
end

cm = cbrewer('seq', 'Blues',10);

H1= clustergram(model.emissions,'Colormap',cm...
    ,'ColumnLabels',markNames...
    ,'Standardize','none'...
    ,'Symmetric', 0);
%     ,'DisplayRange', [0 1]...
adjustPlot
addTitle(H1,'Emissions');

H12 = HeatMap(model.emissions,'Colormap',cm...
    ,'Standardize','none'...
    ,'Symmetric', 0 ...
    ,'ColumnLabels', markNames ...
    ,'RowLabels', 1:param.nB);
%     adjustPlot
addTitle(H12,'Emissions');

if size(model.transitionD,1)>1
    for index = 1:param.nD
        H2 = HeatMap(model.transitionB(param.nB:-1:1,:,index),'Colormap',cm...
            ,'Standardize','none'...
            ,'Symmetric', 0 ...
            ,'ColumnLabels', 1:param.nB ...
            ,'RowLabels', param.nB:-1:1);
        %     adjustPlot
        addTitle(H2,['Nucleosome-Level Transitions in Domain ' num2str(index)]);
    end
else
    H2 = HeatMap(model.transitionB(param.nB:-1:1,:),'Colormap',cm...
        ,'Standardize','none'...
        ,'Symmetric', 0 ...
        ,'ColumnLabels', 1:param.nB ...
        ,'RowLabels', param.nB:-1:1);
    addTitle(H2,'Nucleosome-Level Transitions in Domain 1');
end

if size(model.transitionD,1)>1
    H2 = HeatMap(model.transitionD(param.nD:-1:1,:),'Colormap',cm...
        ,'Standardize','none'...
        ,'Symmetric', 0 ...
        ,'ColumnLabels', 1:param.nD ...
        ,'RowLabels', param.nD:-1:1);
    %     adjustPlot
    addTitle(H2,'Domain Transitions');
end


end

function adjustPlot
% Help from
% http://www.mathworks.com/support/solutions/en/data/1-8E9WQG/index.html?solution=1-8E9WQG

set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');


cgfigidx = strmatch('Clustergram',get(allhnds,'Tag'));
cffighnd = allhnds(cgfigidx);
set(0,'showhiddenHandles','off')
% if length(cffighnd)>1
% warning('More than one clustergram handle found. Using most recent clustergram')
cffighnd = cffighnd(1);
% end

set(cffighnd,'position', [613 64 474 702])
end


%  set(0,'ShowHiddenHandles','on')

%  figureHandle = gcf;

%# make all text in the figure to size 20 and bold
%  set(findall(figureHandle,'type','text'),'fontSize',20,'fontWeight','bold')

%# make all text in the figure to size 20
%  set(findall(gcf,'type','text'),'fontSize',20)
%  set(gcf,'Renderer', 'painters')
