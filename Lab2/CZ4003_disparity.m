% Function to compute disparity map based on two input images and template
% dimensions

% img1: main input image, PI
% img2: secondary input image, PR
% temp_x: template dimensions in x-axis (+- center pixel)
% temp_y: template dimensions in y-axis (+- center pixel)

function [disp_map] = CV4003_disparity(img1, img2, template_x_full, template_y_full)

%% Pre-processing
% Make assertions about input arguments
assert(isnumeric(img1), 'First input image is non-numeric');
assert(isnumeric(img2), 'Second input image is non-numeric');
assert(isnumeric(template_x_full), 'First input dimension is non-numeric');
assert(isnumeric(template_y_full), 'Second input dimension is non-numeric');

% Convert images to double
if ~isa(img1, 'double')
    img1 = double(img1);
end

if ~isa(img2, 'double')
    img2 = double(img2);
end

% Convert template dimensions to +- non-negative integer values
template_x = floor(abs(template_x_full)/2);
template_y = floor(abs(template_y_full)/2);

% Get information about input image sizes
dims_img = size(img1);

% Resize img2 to img1's size if required (not necessary given the two test
% cases we're provided, but doesn't hurt)
if ~isequal(size(img2), dims_img)
    img2 = imresize(img2, dims_img);
end

% Create template for disparity map
disp_map = zeros(dims_img);

% Pad both images with zeros to simplify later operations
img1_pad = zeros(dims_img + [2 * template_x, 2 * template_y]);
img2_pad = zeros(dims_img + [2 * template_x, 2 * template_y]);
img1_pad(template_x + 1 : end - template_x, template_y + 1 : end - template_y) = img1;
img2_pad(template_x + 1 : end - template_x, template_y + 1 : end - template_y) = img2;

%% Computing Disparity Map
% Iterate over all pixels in img1, by looping over rows then columns
for i = 1 : dims_img(1) % Rows
    
    for j = 1 : dims_img(2) % Columns
        
        %% Extracting template from img1
        % Get equivalent coords in the padded image
        padded_i_1 = i + template_y;
        padded_j_1 = j + template_x;
        
        % Extract template
        template_1 = img1_pad(padded_i_1 - template_y : padded_i_1 + template_y, padded_j_1 - template_x : padded_j_1 + template_x);    
        
        %% Computing SSD
        % Set scan line (row)
        padded_i_2 = padded_i_1;
        
        % Set initial SSD value (only changes when i, j changes)
        SSD_base = sum(template_1.^2, 'all');
        
        % Set baseline values
        % 1: SSD value, 2: disparity value
        baseline = [inf, 0];
            
        % Iterate for disparity < 15
        for k = -14 : 14
            
            % Get current column
            padded_j_2 = padded_j_1 + k;
            
            % Check value to prevent overshooting boundaries
            if (padded_j_2 <= template_y) || (padded_j_2 >= dims_img(2) - template_y)
                continue
            end
            
            % Get current template from img2
            template_2 = img2_pad(padded_i_2 - template_y : padded_i_2 + template_y, padded_j_2 - template_x : padded_j_2 + template_x);
            
            % Compute SSD additional terms
            SSD = SSD_base + sum(template_2 .^ 2, 'all') - sum(2 * template_1 .* template_2, 'all');
            
            % Assign new values of SSD and disparity if new SSD is less
            % than current SSD
            if SSD < baseline(1)
                baseline = [SSD, -k]; % k is delta wrt 2 minus 1, when actual disparity should be 1 minus 2
            end
        end
        
        % Assign value to disparity map
        disp_map(i, j) = baseline(2);
    end
end

