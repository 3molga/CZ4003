% Matlab script for CZ4003 Lab 2 Assignment

clear
clc
close all
% cd(''); % This directory is presumably where the code and images are located

%% Edge Detection
% Reading and displaying the image
macritchie_img = imread('.\Images\macritchie.jpg');
macritchie_img_bw = rgb2gray(macritchie_img);
figure, imshow(macritchie_img_bw), title('Grayscale original image');

% Creating Sobel masks and filtering
sobel_horz_mask = [-1 -2 -1; 0 0 0; 1 2 1];
sobel_vert_mask = sobel_horz_mask';
macritchie_vert = conv2(macritchie_img_bw, sobel_vert_mask);
macritchie_horz = conv2(macritchie_img_bw, sobel_horz_mask);
figure, subplot(1, 2, 1), imagesc(macritchie_vert), axis image, axis off, colormap gray;
title('Vertical edges of image');
subplot(1, 2, 2), imagesc(macritchie_horz), axis image, axis off, colormap gray;
title('Horizontal edges of image');

% Combined edge image
macritchie_edge_square = macritchie_vert.^2 + macritchie_horz.^2;
macritchie_edge_square = macritchie_edge_square(3:end-2, 3:end-2); % Remove edges of image that weren't convolved properly
figure, imagesc(macritchie_edge_square), axis image, axis off, colormap gray;
title('Combined edges of image');

% Thresholded edge images
thresholds = [10000, 50000, 100000, 200000];
figure
for i = 1 : length(thresholds)
    macritchie_edge_binary  = macritchie_edge_square > thresholds(i);
    subplot(1, length(thresholds), i), imshow(macritchie_edge_binary), axis image, axis off;
    title(sprintf('Threshold: %01d', thresholds(i)));
    hold on
end

% Canny algorithm
sigmas = [1, 2.3, 3.6, 5];
tls = [0.01, 0.03, 0.06, 0.09];
macritchie_img_bw_d = double(macritchie_img_bw);
macritchie_edge_canny = edge(macritchie_img_bw_d, 'canny', [0.04 0.1], 1);
figure, imagesc(macritchie_edge_canny), axis image, axis off, colormap gray;
title(sprintf('Canny edge image with original parameter values'));

% Varying sigma
figure
for j = 1 : length(sigmas)
    macritchie_edge_canny = edge(macritchie_img_bw_d, 'canny', [0.04 0.1], sigmas(j));
    subplot(1, length(sigmas), j), imagesc(macritchie_edge_canny), axis image, axis off, colormap gray;
    title(sprintf('Sigma: %.1f', sigmas(j)));
    hold on
end

% Varying lower threshold
figure
for k = 1 : length(tls)
    macritchie_edge_canny = edge(macritchie_img_bw_d, 'canny', [tls(k) 0.1], 1);
    subplot(1, length(tls), k), imagesc(macritchie_edge_canny), axis image, axis off, colormap gray;
    title(sprintf('Lower threshold: %.2f', tls(k)));
    hold on
end

%% Line Finding
% Getting image to use
line_img = edge(macritchie_img_bw_d, 'canny', [0.04 0.1], 1);

% Applying Radon transform
[radon_img, xp] = radon(line_img);
figure, subplot(1, 2, 1), imagesc(radon_img), axis image, axis off, colormap gray;
title('Radon transform of edges');

% Comparing versus Hough transform
[hough_img, theta_h, rho_h] = hough(line_img);
subplot(1, 2, 2), imagesc(hough_img), axis image, axis off, colormap gray;
title('Hough transform of edges');

% Getting maximum pixel intensity coordinates
max_pix = max(radon_img, '', 'all');
[radius, theta] = find(radon_img == max_pix);
radius = xp(radius); % Convert back to radial coordinates (for actual distance of radius)

% Deriving normal line equation
[A, B] = pol2cart(theta * pi/180, radius); 
B = -B;
C = A*(A + 179) + B*(B + 145);

% Computing corresponding values for xl and xr
xl = 0;
xr = size(line_img, 2) - 1;
yl = (C - A*xl)/B;
yr = (C - A*xr)/B;

% Displaying overlaid image
figure, imshow(macritchie_img), axis off,
line([xl, xr], [yl, yr], 'LineWidth', 2);
title('Overlaid main line');

%% 3D Stereo
% Read first pair of images and convert to greyscale
corridor_l = rgb2gray(imread('.\Images\corridorl.jpg'));
corridor_r = rgb2gray(imread('.\Images\corridorr.jpg'));

% Display images
figure, subplot(1, 2, 1), imshow(corridor_l), title('Left perspective');
subplot(1, 2, 2), imshow(corridor_r), title('Right Perspective');

% Calculate and display disparity map
corridor_disp = CV4003_disparity(corridor_l, corridor_r, 11, 11);
corridor_disp_ref = imread('.\Images\corridor_disp.jpg');
figure, subplot(1, 2, 1), imagesc(corridor_disp), axis image, axis off, 
title('Computed disparity map');
subplot(1, 2, 2), imagesc(corridor_disp_ref), axis image, axis off, title('Reference disparity map');

% Repeat for triclops image
triclops_l = rgb2gray(imread('.\Images\triclopsi2l.jpg'));
triclops_r = rgb2gray(imread('.\Images\triclopsi2r.jpg'));
figure, subplot(1, 2, 1), imshow(triclops_l), title('Left perspective');
subplot(1, 2, 2), imshow(triclops_r), title('Right perspective');

triclops_disp = CV4003_disparity(triclops_l, triclops_r, 11, 11);
triclops_disp_ref = imread('.\Images\triclopsid.jpg');
figure, subplot(1, 2, 1), imagesc(triclops_disp), axis image, axis off, title('Computed disparity map');
subplot(1, 2, 2), imagesc(triclops_disp_ref), axis image, axis off, title('Reference disparity map');
