% Matlab script for CZ4003 Lab 1 Assignment

clear
clc
close all
% cd(''); % This directory is presumably where the code and images are located

%% Part 1: Contrast stretching
% Read and display original image
Pc = imread('.\images\mrt-train.jpg');
whos Pc
P = rgb2gray(Pc);
figure, subplot(1, 2, 1), imshow(P), title('Original MRT Train Image')

% Getting min and max intensities of image
P_min = uint8(min(P(:)));
P_max = uint8(max(P(:)));

% Contrast stretching
P2 = P - P_min;
P2 = P2 .* (255/double(P_max - P_min));

% Min and max intensities of contrast-stretched image, and
P2_min = uint8(min(P2(:)));
P2_max = uint8(max(P2(:)));
subplot(1, 2, 2), imshow(P2), title('Contrast-stretched MRT Train Image');

%% Part 2: Histogram Equalization
% Displaying histograms with 10 and 256 bins
figure, sgtitle('Original Image''s Histograms') 
subplot(1, 2, 1), imhist(P, 10), title('10-bin Histogram');
subplot(1, 2, 2), imhist(P, 256), title('256-bin Histogram');

% Histogram equalization - first round
P3 = histeq(P, 255);
figure, sgtitle('Equalized Image''s Histograms') 
subplot(1, 2, 1), imhist(P3, 10), title('10-bin Histogram');
subplot(1, 2, 2), imhist(P3, 256), title('256-bin Histogram');

% Histogram equalization - second round
P4 = histeq(P3, 255);
figure, sgtitle('2nd Round Equalized Image''s Histograms') 
subplot(1, 2, 1), imhist(P4, 10), title('10-bin Histogram');
subplot(1, 2, 2), imhist(P4, 256), title('256-bin Histogram');

%% Part 3: Linear Spatial Filtering
% Creating Gaussian filters
gauss1 = fspecial('gaussian', 5, 1);
gauss2 = fspecial('gaussian', 5, 2);

% Viewing the Gaussian filters
figure, sgtitle('Gaussian filter 3D view')
subplot(1, 2, 1), mesh(gauss1), title('\sigma = 1');
subplot(1, 2, 2), mesh(gauss2), title('\sigma = 2');

% Reading, displaying and filtering the image with additive Gaussian noise
gaussNoiseImg = imread('.\images\lib-gn.jpg');
gaussNoiseImgConv1 = conv2(gaussNoiseImg, gauss1); % First filter
gaussNoiseImgConv2 = conv2(gaussNoiseImg, gauss2); % Second filter
figure, subplot(1, 3, 1), imshow(gaussNoiseImg), title('Original Image')
subplot(1, 3, 2), imagesc(gaussNoiseImgConv1); axis image, axis off, 
title('Filtered with first Gaussian filter')
subplot(1, 3, 3), imagesc(gaussNoiseImgConv2); axis image, axis off, 
title('Filtered with second Gaussian filter'), colormap gray;

% Reading, displaying and filtering the image with additive speckle noise
speckleNoiseImg = imread('.\images\lib-sp.jpg');
speckleNoiseImgConv1 = conv2(speckleNoiseImg, gauss1);
speckleNoiseImgConv2 = conv2(speckleNoiseImg, gauss2);
figure, subplot(1, 3, 1), imshow(speckleNoiseImg), title('Original Image');
subplot(1, 3, 2), imagesc(speckleNoiseImgConv1); axis image, axis off, 
title('Filtered with first Gaussian filter')
subplot(1, 3, 3), imagesc(speckleNoiseImgConv2); axis image, axis off, 
title('Filtered with second Gaussian filter'), colormap gray;

%% Part 4: Median Filtering
% Repeating steps for image with additive Gaussian noise
gaussNoiseImgMed1 = medfilt2(gaussNoiseImg, [3 3]);
gaussNoiseImgMed2 = medfilt2(gaussNoiseImg, [5 5]);
figure, subplot(1, 3, 1), imshow(gaussNoiseImg), title('Original Image');
subplot(1, 3, 2), imagesc(gaussNoiseImgMed1); axis image, axis off
title('Filtered with 3x3 median filter')
subplot(1, 3, 3), imagesc(gaussNoiseImgMed2); axis image, axis off
title('Filtered with 5x5 median filter'), colormap gray;

% Repeating steps for image with additive speckle noise
speckleNoiseImgMed1 = medfilt2(speckleNoiseImg, [3 3]);
speckleNoiseImgMed2 = medfilt2(speckleNoiseImg, [5 5]);
figure, subplot(1, 3, 1), imshow(speckleNoiseImg), title('Original Image');
subplot(1, 3, 2), imagesc(speckleNoiseImgMed1); axis image, axis off
title('Filtered with 3x3 median filter')
subplot(1, 3, 3), imagesc(speckleNoiseImgMed2); axis image, axis off
title('Filtered with 5x5 median filter'), colormap gray;

%% Part 5: Supressing Noise Interference Problems
% Reading and displaying the image
pck_img = imread('.\images\pck-int.jpg');
figure, imshow(pck_img), title('Original image');

% Fourier transform
pck_F = fft2(pck_img);
pck_S = abs(pck_F);
figure, subplot(1, 2, 1), imagesc(fftshift(pck_S.^0.1)); colormap('default'), axis image 
title('Shifted power spectrum of original Fourier transform');
subplot(1, 2, 2), imagesc(pck_S.^0.1); colormap('default'), axis image
title('Original power spectrum of original Fourier transform');
% [pck_F_x, pck_F_y] = ginput(2); % Rough values below
pck_F_x =  [9, 249];
pck_F_y = [240, 16];

% Setting neighbourhood elements to zero and recomputing
pck_F2 = pck_F;
for i = 1 : 2
    cur_points = [round(pck_F_x(i)), round(pck_F_y(i))];
    pck_F2(cur_points(2)-2 : 1 : cur_points(2)+2, cur_points(1)-2 : 1 : cur_points(1)+2) = 0;
end
pck_S2 = abs(pck_F2);
figure, subplot(1, 2, 1), imagesc(fftshift(pck_S2.^0.1)); colormap('default'), axis image 
title('Shifted power spectrum of altered Fourier transform');
subplot(1, 2, 2), imagesc(pck_S2.^0.1); colormap('default'), axis image
title('Original power spectrum of altered Fourier transform');

% Inverse fourier transform
pck_img2 = real(uint8(ifft2(pck_F2)));
figure, imshow(pck_img2), title('Result of altered Fourier transform');

% Trying to improve the result
% How about attenuating additional rows/columns of the frequency spectrum
% Expanding the 5x5 neighbourhood to 7x7
% And changing from zeroing to attenutating to the 25th percentile
pck_F3 = pck_F;
pck_25 = prctile(pck_F3, 25, 'all');
for i = 1 : 2
    cur_points = [round(pck_F_x(i)), round(pck_F_y(i))];
    pck_F3(cur_points(2)-3 : 1 : cur_points(2)+3, cur_points(1)-3 : 1 : cur_points(1)+3) = pck_25;
    pck_F3(cur_points(2), :, :) = pck_25; % Set entire row to 25th percentile
    pck_F3(:, cur_points(1), :) = pck_25; % Set entire column to 25th percentile
end
pck_img3 = real(uint8(ifft2(pck_F3)));
figure, subplot(1, 2, 1), imshow(pck_img3), title('Attempt at improving result further');
subplot(1, 2, 2), imagesc(fftshift(abs(pck_F3).^0.1)), axis image
title('Power spectrum of further altered Fourier transform');

%% Part 5.1 Other Noise Interference Problem
%Reading other primate image
primate_img = imread('.\images\primate-caged.jpg');
primate_img = rgb2gray(primate_img);
figure, subplot(1, 2, 1), imshow(primate_img), title('Original Image')

% Fourier transform
primate_F = fft2(primate_img);
primate_S = abs(primate_F);
subplot(1, 2, 2), imagesc(fftshift(primate_S.^0.1)), axis image 
title('Power spectrum of original Fourier transform');

% Instead of specifying points, what about specifying a ROI to attenuate
polygon = drawpolygon;
mask = createMask(polygon);

% Create 11x11 neighbourhood as exception to attenuation centred on [129, 129]
%([0, 0] on non-shifted Fourier transform - DC/near-DC components)
mask(124:134, 124:134) = 0;

% Attenuate power spectrum of Fourier transform
primate_F2 = fftshift(primate_F);
primate_F2(mask) = primate_F2(mask) .^ 0.7; % Arbitrary parameter
primate_img2 = real(uint8(ifft2(ifftshift(primate_F2))));
figure, subplot(1, 3, 1), imshow(primate_img2), title('Attempt at freeing the primate');
subplot(1, 3, 2), imagesc(abs((primate_F2)).^0.01)
title('Power spectrum of altered Fourier transform'), axis image % Show attenuated Fourier transform

% Attempting to improve the result by applying Gaussian filter from previous section
primate_img3 = imfilter(primate_img2, gauss1, 'conv', 'replicate');
subplot(1, 3, 3), imshow(primate_img3), title('Attempt at improving the previous image');
% It's... different, at least

%% Part 6: Undoing Perspective Distortion
% Reading, displaying and getting coords of book corners
book_img = imread('.\images\book.jpg');
figure, subplot(1, 2, 1), imshow(book_img), title('Original image')
[book_ori_x, book_ori_y] = ginput(4); % Clockwise from bottom left corner

% Setting coords of transformed book corners
book_transform_x = [0; 0; 210; 210]; % Assuming 210 x 297 px
book_transform_y = [297; 0; 0; 297];

% Setting up transformed point vector
book_transform_vec = zeros(8, 1); % Preallocating
vec_point = 1;
for n_point = 1 : 4 % Loop over each point in transformed image (which then has 2 coordinate values)
    % x-coords
    book_transform_vec(vec_point) = book_transform_x(n_point);
    vec_point = vec_point + 1;
    
    % y-coords
    book_transform_vec(vec_point) = book_transform_y(n_point);
    vec_point = vec_point + 1;
end

% Setting up projection matrix thing
% Each row is a row vector with 8 entries
for n_point = 1 : 4
    % Two rows for each , just fill them in directly I suppose
    for row = 1 : 2
        temp_row = zeros(1, 8); % Set up temporary row vector
        if mod(row, 2) % 1st row
            temp_row(1:2) = [book_ori_x(n_point), book_ori_y(n_point)];
            temp_row(3) = 1;
            temp_row(7) = -(book_transform_x(n_point)).*temp_row(1);
            temp_row(8) = -(book_transform_x(n_point)).*temp_row(2);
        else % 2nd row
            temp_row(4:5) = [book_ori_x(n_point), book_ori_y(n_point)];
            temp_row(6) = 1;
            temp_row(7) = -(book_transform_y(n_point)).*temp_row(4);
            temp_row(8) = -(book_transform_y(n_point)).*temp_row(5);
        end
        
        % Add temporary row to projection matrix
        if n_point == 1 && row == 1
            proj_mat = temp_row;
        else
            proj_mat = cat(1, proj_mat, temp_row);
        end
    end
end

% Calculations
transform_params = proj_mat \ book_transform_vec;
transform_params_mat = reshape([transform_params; 1], 3, 3)';

% Verification
transformed_points = transform_params_mat * [book_ori_x'; book_ori_y'; ones(1, 4)];
transformed_points = transformed_points ./ (ones(3, 1) * transformed_points(3, :));
transformed_points = round(transformed_points);

% Warping
T_warp = maketform('projective', transform_params_mat');
book_img_T = imtransform(book_img, T_warp, 'XData', [0 210], 'YData', [0 297]);
subplot(1, 2, 2), imagesc(book_img_T), title('Transformed image'), axis image, axis off

%% Attempting to identify the pink area
% Getting mask of all pink-ish areas in YCbCr colorspace
book_ycbcr = rgb2ycbcr(book_img_T);
mask = book_ycbcr(:, :, 3) > 131; % Semi-arbitrary but effective threshold
book_mask1 = labeloverlay(book_img_T, mask, 'Transparency', 0.4);
figure, subplot(2, 2, 1), imshow(book_mask1), axis image, title('Original mask')

% Mask is kind of noisy/dirty, use morphological opening to clean it up somewhat
strel_open = strel("disk", 3);
mask_opened = imopen(mask, strel_open); % Morphological opening
book_mask2 = labeloverlay(book_img_T, mask_opened, 'Transparency', 0.4);
subplot(2, 2, 2), imshow(book_mask2), axis image
title('Mask after morphological opening')

% Filter out all areas in the mask except for the one with the largest area
mask_regions = regionprops(mask_opened, 'Area', 'PixelIdxList'); % Regionprops call
[~, mask_region_maxID] = max([mask_regions.Area]); % Get index of region with largest area
mask_regions = mask_regions(mask_region_maxID); % Remove all other regions from regionprops struct
mask_filtered = false(size(mask_opened)); % Recreate mask
mask_filtered([mask_regions.PixelIdxList]) = 1;
book_mask3 = labeloverlay(book_img_T, mask_filtered, 'Transparency', 0.4);
subplot(2, 2, 3), imshow(book_mask3), axis image
title('Mask after filtering for largest continuous area')

% Do morphological closing to fill in the part of the mask that is "dented"
strel_close = strel("disk", 20);
mask_closed = imclose(mask_filtered, strel_close); % Morphological closing
book_mask4 = labeloverlay(book_img_T, mask_closed, 'Transparency', 0.4);
subplot(2, 2, 4), imshow(book_mask4), axis image
title('Mask after morphological closing');