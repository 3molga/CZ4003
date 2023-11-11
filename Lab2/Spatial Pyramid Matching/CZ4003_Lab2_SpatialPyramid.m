% Example script to test feature extraction methods from Lazebnik, Schmid
% and Ponce's 2006 paper "Beyond Bags of Features: Spatial Pyramid Matching
% for Recognizing Natural Scene Categories"

% Due to time and computational constraints, only the functions required to extract the
% features are included, and the models specified in the Bag-of-Words vs
% Spatial Pyramid Matching example aren't trained/tested

% Also, currently assumes that all input images are of a similar size

clear
clc
% addpath(genpath('')); % This directory is presumably where the code and images are located

%% Part 1: Loading A Test Image
% Entire Caltech-101 dataset was downloaded, but only one test image is
% being used here (image name specified from "wild-cat" category)
im_test = imread('image_0010.jpg');

% Part 5 of the paper specifies that all experiments were conducted in
% grayscale
im_test = rgb2gray(im_test);

%% Part 2: Extracting Features
% Two types of features are specified to be extracted in the experiments
% "Weak" features: points with gradient magnitude in a direction exceeding
% a threshold
% "Strong" features: SIFT descriptors of 16x16 patches with spacings of 8
% pixels
% I will attempt to extract both (though no guarantees that it even works
% as intended)

weakFeatures = extractWeakFeatures(im_test);
strongFeatures = extractStrongFeatures(im_test);

%% Part 3: Theoretical code to train SVMs
%% Partition train and test set
% First get list of all images
% Replace with path to Caltech-101 dataset that has been downloaded
% img_list = recursiveListFiles('', '.jpg');

% Then identify class, which is specified in the folder name
class_list = arrayfun(@returnFolder, img_list, 'UniformOutput', false);
class_unique = unique(class_list); 

% Then split into train and test sets
% Part 5.2 of the paper specifies that for each class, there are 30
% training images and 50 test images
for i = 1 : numel(class_unique)
    
    % Get current class
    cur_class = class_unique{i};
    
    % Match against class_list to get IDs of images for that class
    cur_IDs = find(ismember(class_list, cur_class));
    
    % Extract 80 random IDs
    % Some classes are too small for even 80 images total, so use a min
    % statement instead
    randIDs = randperm(numel(cur_IDs), 1, min(80, numel(cur_IDs)));
    
    % Split into train and test
    trainIDsCur = cur_IDs(randIDs(1:30));
    testIDsCur = cur_IDs(randIDs(31:end));
    
    % Add to list
    if i == 1
        trainIDs = trainIDsCur;
        testIDs = testIDsCur;
    else
        trainIDs = cat(1, trainIDs, trainIDsCur);
        testIDs = cat(1, testIDs, testIDsCur);
    end
end

% Extract corresponding labels and image paths
train_labels = class_list(trainIDs);
test_labels = class_list(testIDs);
train_imgs = img_list(trainIDs);
test_imgs = img_list(testIDs);

%% Extract Features
% Not actually used due to aforementioned limitations
% Iterate over train sets and extract features
for train_i = 1 : numel(trainIDs)
    
    % Load image
    im = imread(train_imgs(train_i));
    
    % Make sure image is greyscale
    if numel(size(im)) > 2
        im = rgb2gray(im);
    end
    
    % Extract weak and strong features
    weakFeature = extractWeakFeatures(im);
    strongFeature = extractStrongFeatures(im);
    
    % Extract random patches from strongFeature
    % As specified in section 4 of the paper, but no details on the
    % granularity/number of patches were specified, so I'll assume that 10%
    % of the patches are chosen
    randpatchIDs = randi(numel(strongFeature), floor(0.1*numel(strongFeature)), 1);
    patches = strongFeature(randpatchIDs);
    
    % Concatenate features together
    if train_i == 1
        totalWeak = weakFeature;
        totalStrong = patches;
    else
        totalWeak = cat(1, totalWeak, weakFeature);
        totalStrong = cat(1, totalStrong, patches);
    end
end