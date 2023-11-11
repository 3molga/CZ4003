% Function to extract "strong" features from an input image as per the 2006 paper

% From my understanding, "strong" features are SIFT descriptors of 16x16
% patches with 8 pixel spacings in-between
% The paper doesn't specify how they deal with inconsistent image sizes, so
% I'll just assume we step until it is no longer possible to step by
% another 8 pixels

% Either the original text doesn't describe the features in enough detail,
% or I'm just not smart enough for this, so there's a very high chance that
% the method I'm implementing is completely wrong (how do they get M =
% 200/400 from SIFT detectors that should be 128-long?)

function [features] = extractStrongFeatures(inputImg)

% Pre-processing
% Start by calculating number of steps needed
imSize = size(inputImg);
stepCounts = floor(imSize./8) - [1, 1];

% Pre-allocte cell array for SIFT descriptors
features = cell(stepCounts);

% Step along image grid and extract corresponding SIFT descriptors
for i = 1 : stepCounts(1)
    
    for j = 1 : stepCounts(2)
        
        % Get current patch
        i_center = i * 8;
        j_center = j * 8;
        imPatch = inputImg(i_center - 7 : i_center + 8, j_center - 7 : j_center + 8);
        
        % Get SIFT descriptor
        SIFTdescriptor = extractSIFTdescriptor(imPatch);
        
        % Concatenate as cells
        features{i, j} = SIFTdescriptor;
        
    end
    
end

end

