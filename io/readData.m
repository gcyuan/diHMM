function emissions = readData(cellType, chr, param)
% Read binarized data

dataBinaryFile = fullfile(param.fullBinDataFolder.(cellType), [cellType '_' chr '_binary.txt']);
dataBin2decFile = fullfile(param.fullBin2decDataFolder.(cellType), [cellType '_' chr '_bin2dec.txt']);

% We read binary data
fid = fopen(dataBinaryFile);
% First line should contain cellLine and chr
s = fgetl(fid);
I = strfind(s, char(9));

cellTypeFromFile = s(1:I(1)-1);
chrFromFile = s(I(1)+1:end);
if ~strcmp(cellType, cellTypeFromFile) || ~strcmp(chr, chrFromFile)
    error('Reading the wrong file it seems!')
end

% Second line should contain the list of histone marks
s = fgetl(fid);
I = strfind(s, char(9));

nMarks = length(I)+1;
markNames = cell(1, nMarks);
if nMarks > 1
    markNames{1} = s(1:I(1)-1);
    for index = 1:nMarks-2
        markNames{index+1} = s(I(index)+1:I(index+1)-1);
    end
    markNames{nMarks} = s(I(nMarks-1)+1:end);
else
    markNames{1} = s;
end

data = textscan(fid, repmat('%d',1,nMarks), 'delimiter', char(9));
fclose(fid);
% numData = [logical(cell2mat(data)) ~logical(cell2mat(data))]; % We save space, just in case
numData = logical(cell2mat(data)); % We save space, just in case
% numData = uint8(cell2mat(data)); % We save space, just in case

% We crop the data to fit into multiples of our domain size
numData = numData(1:end-mod(size(numData,1),param.nDSize),:);

% We read bin2dec data, if ~exist, create first
if exist(dataBin2decFile,'file')
    fid = fopen(dataBin2decFile);
    % First line should contain cellLine and chr
    s = fgetl(fid);
    I = strfind(s, char(9));
    
    cellTypeFromFile = s(1:I(1)-1);
    chrFromFile = s(I(1)+1:end);
    if ~strcmp(cellType, cellTypeFromFile) || ~strcmp(chr, chrFromFile)
        error('Reading the wrong file it seems!')
    end
    
    % Second line should contain the list of histone marks
    s = fgetl(fid);
    I = strfind(s, char(9));
    
    nMarks = length(I)+1;
    markNames = cell(1, nMarks);
    if nMarks >1
        markNames{1} = s(1:I(1)-1);
        for index = 1:nMarks-2
            markNames{index+1} = s(I(index)+1:I(index+1)-1);
        end
        markNames{nMarks} = s(I(nMarks-1)+1:end);
    else
        markNames{1} = s;
    end
    
    bin2decData = cell2mat(textscan(fid, '%f'));
    fclose(fid);
else
    bin2decData = bin2dec(num2str(cell2mat(data)));
    bin2decData = bin2decData(1:size(numData,1));
    if ~exist(fileparts(dataBin2decFile), 'dir') % Is bin2dec directory there?
        mkdir(fileparts(dataBin2decFile))
    end
    
    fs = fopen(dataBin2decFile, 'w+');
    header1 = {cellType chr};
    
    header2 = markNames;
    
    format2 = [repmat('%s\t',1,nMarks-1), '%s\n'];
    fprintf(fs, '%s\t%s\n', header1{:});
    fprintf(fs, format2, header2{:});
    fprintf(fs, '%d\n', bin2decData);
    fclose(fs);
end
    
emissions.markNames = markNames;
emissions.binData = numData;
emissions.bin2decData = bin2decData;

% emissions.binTopBoundary = param.binSize*(1:size(numData,1));