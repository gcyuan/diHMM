function chrList = getChrList(genome)

if ~exist('genome', 'var')
    genome = 'hg19';
end

switch genome
    case {'hg19', 'hg18', 'hg17'}
        chrList = cell(23,1);
        for index = 1:length(chrList)-1
            chrList{index} = ['chr' num2str(index)];
        end
        chrList{23} = 'chrX';
    case 'mm9'
        chrList = cell(20,1);
        for index = 1:length(chrList)-1
            chrList{index} = ['chr' num2str(index)];
        end
        chrList{20} = 'chrX';
end
