
%% Read the image
img = imread("~/Desktop/TFM/bsc83671/GTEx_v8/Laura/00.Data/Fem_tissues_images/Vagina/GTEX-ZAJG-2925.svs",1);
disp(size(img));
img = img(:, :, 1:3);
%%
imshow(img);

%% Initialize a mask with zeros (same size as the image) - for 2D grayscale or RGB images
mask = false(size(img, 1), size(img, 2));  % Binary mask, same size as the image

%% Loop to select multiple regions
while true

    % Select a region using freehand drawing
    h = drawfreehand;
    
    % Create a binary mask from the selected region

    region_mask = createMask(h);    
    % Combine the new mask with the overall mask using logical OR
    mask = mask | region_mask;
    
    % Ask the user if they want to select another region
    answer = questdlg('Do you want to select another region?', ...
                      'Continue Selection', 'Yes', 'No', 'Yes');
    
    % If the user selects "No", break the loop
    if strcmp(answer, 'No')
        break;
    end
end

%% Display the combined 
figure;
imshow(mask);
title('Combined Binary Mask');
%%
% Resize the mask to a smaller size, e.g., 50% of the original size
scaled_mask = imresize(mask, 0.5); % You can adjust the scaling factor as needed

% Save the resized mask as a JPEG image
imwrite(scaled_mask, '/home/lventura/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/Segmentation/FallopianTube/vessels/GTEX-NPJ7-1225.jpg');

%% Overlay the mask on the original image
imwrite(mask,'~/Desktop/TFM/bsc83671/GTEx_v8/Laura/12.Tissues_substructures/Segmentation/Uterus/endometrium/GTEX-1ICG6-1425.jpg');


