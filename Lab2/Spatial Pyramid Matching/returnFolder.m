% Simple function to return the folder from the full path of a file

function folder = returnFolder(path)

% Split based on file seps
split = strsplit(path, filesep);

% Folder is one level before last entry (which is file name)
folder = split{end - 1};
end
