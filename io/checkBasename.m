function checkBasename(baseName)
% We check if dir exists, if not we create it

[pathstr, ~, ~] = fileparts(baseName);

if ~exist(pathstr,'dir')
    mkdir(pathstr)
end
end