function chromSizes = getChromSizes(file)
% get chromosome sizes for reference genome
fid = fopen(file, 'r');
data = textscan(fid, '%s %d');    % It has two columns
fclose(fid);

for index = 1:length(data{1})
    chromSizes.(data{1}{index}) = data{2}(index);
end
