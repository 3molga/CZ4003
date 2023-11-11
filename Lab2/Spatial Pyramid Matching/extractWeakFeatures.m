% Function to extract "weak" features from an input image as per the 2006 paper

% From my understanding, "weak" features refer to points with gradient
% magnitude in one of 8 directions (likely 45deg angles) exceeding a
% specified threshold (the details of which aren't quite specified)
% As such, I'll just set a "global" threshold of being above the 95th
% percentile over all gradient magnitudes in all directions

% Returns 1 cell array containing 16 (channels) points as cells, 
% each cell containing a matrix of n (number of points) x 2 (row and col coords)

function [features] = extractWeakFeatures(inputImg)

% Pre-process image if necessary
imSize = size(inputImg);

if length(imSize) > 2 % Image isn't grayscale
    inputImg = rgb2gray(inputImg);
end

k = 1;
features = cell(16, 1);

% Loop over two image scales
% First scale: original image
% Second scale: image downsampled by factor of 2
for scale = 1 : 2
    if scale == 2
        inputImg = inputImg(1:2:end, 1:2:end);
    end
    
    % Preallocate 3D matrix to store gradient magnitudes
    grad_mags = zeros([ceil(imSize ./ scale), 8]);
    
    % Loop over all directions and compute histograms
    for i = 1 : 8
        
        % There doesn't seem to be an elegant way to compute gradients in
        % pre-specified directions in Matlab? So I'll just hack something
        % together that approximates the Kirsch operator but with Sobel kernels
        % https://en.wikipedia.org/wiki/Kirsch_operator
        % Basically rotates from facing up by 45deg clockwise
        kernel = returnKernelRotated(i);
        
        % Apply kernel
        imFilt = imfilter(inputImg, kernel, 'replicate');
        
        % Add to grad_mags
        grad_mags(:, :, i) = imFilt;
    end
    
    % Compute 95th percentile cut-off
    perc95 = prctile(grad_mags, 95, 'all');
    
    % Get I, J coordinates from each kernel  
    for j = 1 : 8
        [features_row, features_col] = find(grad_mags(:, :, j) >= perc95);
        
        % Assign to cell 
        features{k} = [features_row, features_col];
        
        % Increment k
        k = k + 1;
    end
end

end

