% Function to extract SIFT detector from 16x16 patch per the 2006 paper
% Follows the general method described in http://www.scholarpedia.org/article/Scale_Invariant_Feature_Transform#Image_descriptor

function [SIFThistogram] = extractSIFTdescriptor(inputImg)

% Make assertion about input image size
assert(isequal(size(inputImg), [16 16]), 'Input image patch should be 16x16!');

% Compute patch gradient and magnitudes
[gradMag, gradDir] = imgradient(inputImg);

% Weigth gradient magnitudes by Gaussian window function (standard sigma = 2)
gradMag = gradMag .* fspecial('gaussian', [16 16], 2);

% Base direction bin threshold
baseBinDir = [-45/2, 45/2]; % Corresponds to gradient direction vertically up

% Compute SIFT detector
for i_start = 0 : 3 % Iterate over rows
    
    for j_start = 0 : 3 % Iterate over columns
        
        % Preallocate current 4x4 patch histogram
        SIFTCur = zeros(1, 8);
        
        % Get top left corner of grid
        topleft = [1 + 4 * i_start, 1 + 4 * j_start];
        
        % Extract grid from gradMag and gradDir
        gradMagCur = gradMag(topleft(1) : topleft(1) + 3, topleft(2) : topleft(2) + 3);
        gradDirCur = gradDir(topleft(1) : topleft(1) + 3, topleft(2) : topleft(2) + 3);
        
        for dir = 0 : 7 % Iterate over bins of directions
            
            % Get current direction bin thresholds
            dirCur = baseBinDir + dir * 45;
            
            % Get magnitude entries whose direction falls within the
            % thresholds
            gradMagCurBinned = gradMagCur((dirCur(1) < gradDirCur) & (gradDirCur < dirCur(2)));
            
            % Add special exception cases where the direction happens to
            % fall on the threshold itself
            gradMagCurBinnedSpecial = gradMagCur(gradDirCur == dirCur(1) | gradDirCur == dirCur(2));
            
            % Get sum of magnitudes
            gradMagSum = sum(gradMagCurBinned, 'all');
            
            % If any special cases are found, increment the sum by half
            % their magnitude values (since they'll also be included in the
            % bin of the prev/next direction)
            if any(gradMagCurBinnedSpecial)
                gradMagSum = gradMagSum + sum(gradMagCurBinnedSpecial, 'all')/2;
            end
            
            % Allocate to SIFTCur
            SIFTCur(dir + 1) = gradMagSum;
        end
        
        % Allocate to overall histogram
        if (i_start == 0) && (j_start == 0)
            SIFThistogram = SIFTCur;
        else
            SIFThistogram = cat(2, SIFThistogram, SIFTCur);
        end
    end

end

end
