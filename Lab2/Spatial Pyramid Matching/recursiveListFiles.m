% Function to recursively list files that match an expression in a certain
% directory

function [fileList] = recursiveListFiles(directory, pattern)

% Make assertions
assert(ischar(directory), 'Directory should be specified as char');
assert(ischar(pattern), 'Pattern should be specified as char');

% Move to directory
oldDir = cd(directory);

% List all entries in directory
tempList = ls;

% Create temporary file list
fileList = strings(0);

% Go over all entries (skipping first 2 for . and ..) and recursively
% search for files matching the pattern
for i = 3 : size(tempList, 1)
    
    % Create path of current entry
    currentEntry = sprintf('%s%s%s', directory, filesep, tempList(i, :));
    currentEntry = strrep(currentEntry, ' ', '');
    
    % If entry does not have a period, it's a folder
    if isfolder(currentEntry)
        newFileList = recursiveListFiles(currentEntry, pattern);
        fileList = cat(1, fileList, newFileList);
    end
    
    % If entry has a period, it's a file
    if isfile(currentEntry)
        if contains(currentEntry, pattern) % On R2019a, would use matches instead if on R2019b+
            fileList = cat(1, fileList, currentEntry);
        end
    end            
end

% Return to previous directory
cd(oldDir);

end

