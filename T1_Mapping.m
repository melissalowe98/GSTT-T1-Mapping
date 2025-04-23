clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD DICOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Select MOLLI Folder set')
dirName = uigetdir(); 
options = struct('recursive', true, 'verbose', true, 'loadCache', false);
[partitions, meta] = readDicomSeries(dirName, options);

% Dynamically calculate the number of TIs
nbti = length(partitions);
fprintf('Number of inversion times (nbti) calculated: %d\n', nbti);

 % Return values:
%   imagePartitions: Array of structs containing all partitions found
%   metaFilenames: Cell array of dicom filenames that contain no images

[image1, info1] = readDicomSeriesImage(dirName, partitions(1));

nbrow = size(image1,1);
nbcol = size(image1,2);
nbslice = size(image1, 3);
nbvols = length(partitions);
nbvoxels = nbrow*nbcol*nbslice;
nbseries = length(partitions);

tinv_acq = zeros(nbti,1);
ttrig_acq = zeros(nbti,1);

for i=1:nbti
    [image, info] = readDicomSeriesImage(dirName, partitions(i));
    metadata(i) = info{1,1}; % to get the dicom headers for every file (TI)
    
    tinv_acq(i)=metadata(i).InversionTime;  % builds a vector of all of the TI's
end

% % need to reorder the TI's in tinv(i), 11 TI's in MOLLI
[tinv,new_order] = sort(tinv_acq);

for k = 1:nbseries
    [image, info] = readDicomSeriesImage(dirName, partitions(k));
	dataTmp = image;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice 
		dataTmp2(:,:,ss,k) = dataTmp(:,:,ss); 
    end
end 
for j = 1:nbseries
    ordernum=new_order(j)
    data(:,:,:,j)= dataTmp2(:,:,:,ordernum);
end

size(data)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Masking to 4D Data (rows x cols x slices x inversion_times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the mask on the first slice and first inversion time
currSlice = squeeze(data(:, :, 1, 1)); % Extract first slice and first TI for mask selection
imshow(currSlice, []);
title('Select ROI (applied to all slices and TIs)');

% Use ROI selection tool (you can change this as needed)
h = imellipse; % Or use roipoly/imfreehand/imrect
wait(h);
BW = createMask(h); % Create binary mask from selected ROI

% Expand the mask to cover all slices and inversion times
roiMask4D = repmat(BW, [1, 1, nbslice, nbti]); % This creates a 4D mask

% Apply the 4D mask to the data
maskedData = data .* roiMask4D;

% Optionally: Set masked (background) data to NaN or any other value
%maskedData(maskedData == 0) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Registration Process for 4D Data (rows x cols x slices x
% inversion_times) - initialise registeredSeries/refIndex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize an empty matrix for the registered images
registeredSeries = zeros(size(maskedData));

% Reference image for registration (use the middle TI as the reference)
%refIndex = round(nbti / 2);
refIndex = round(5); %chose best slice - pick slice with best contrast (ADD CHOOSING STEP?), 6 actually seems to give best results?
refImage = maskedData(:, :, :, refIndex); % Take the entire slice and all slices for the reference

% Copy the reference image into the registered series
registeredSeries(:, :, :, refIndex) = refImage;

% Set up the optimizer and metric for registration
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius / 3;
optimizer.MaximumIterations = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Registration Method: sequential registration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORWARD REGISTRATION (refIndex → nbti)
for movingIndex = refIndex + 1 : nbti
    prevIndex = movingIndex - 1; % Use last registered image as reference

    movingImage = maskedData(:, :, :, movingIndex);
    referenceImage = registeredSeries(:, :, :, prevIndex);

    % Skip if the moving image is empty
    if all(movingImage(:) == 0) || all(referenceImage(:) == 0)
        disp(['Skipping empty slice at inversion time ', num2str(movingIndex)]);
        continue;
    end

    % Perform registration
    registeredSeries(:, :, :, movingIndex) = imregister(movingImage, referenceImage, 'rigid', optimizer, metric);
end

%% BACKWARD REGISTRATION (refIndex → 1)
for movingIndex = refIndex - 1 : -1 : 1
    prevIndex = movingIndex + 1; % Use last registered image as reference

    movingImage = maskedData(:, :, :, movingIndex);
    referenceImage = registeredSeries(:, :, :, prevIndex);

    % Skip if the moving image is empty
    if all(movingImage(:) == 0) || all(referenceImage(:) == 0)
        disp(['Skipping empty slice at inversion time ', num2str(movingIndex)]);
        continue;
    end

    % Perform registration
    registeredSeries(:, :, :, movingIndex) = imregister(movingImage, referenceImage, 'similarity', optimizer, metric);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % T1 FITTING
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

registered_data = registeredSeries;

nbrow2 = size(registered_data,1);
nbcol2 = size(registered_data,2);
nbvoxels2 = nbrow2*nbcol2*nbslice;

% create a mask to speed up calc, thresholds data 300 in final volume
mask=registered_data(:,:,nbti); %shows data in the final TI
mask(le(mask,100))=1; % changed from 300 to 100 2020 July EB
mask(ge(mask,100))=1; % not sure the point in having le and ge the same

% initialise matrices 
t1vec = zeros(1,nbvoxels2,'single'); %Maybe need to reduce this to be new data size
t1vec_star = zeros(1,nbvoxels2,'single');
residuals = zeros(1, nbvoxels2);  % One residual value per voxel

slope = zeros(nbti,nbvoxels2); % there are 11 TI's in the MOLLI sequence

%% Calculate T1

indechs = 1;

% create a 2D array with TI's as the 2nd dimension
for z=1:nbslice
    for y=1:nbcol2
        for x=1:nbrow2  
        if (mask(x,y,z)==1)
            slope(:,indechs) = registered_data(x,y,z,:);  %maybe we do want data to be 4D
        end
        indechs = indechs + 1;
        end
    end 
end

fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'Upper',[6000,12000,5000],'StartPoint',[1000,2000,1000]);

molli = fittype('abs(Axy - Bxy * exp(-tinv/tonestar))','dependent',{'y'},...
    'independent',{'tinv'},'coefficients',{'Axy','Bxy','tonestar'},'options',fo);

for i=1:nbvoxels2
    recover = slope(:,i);
    i
    if recover(1)~=0
                f = fit(tinv,recover,molli); 
                coeffvals = coeffvalues(f);
                Tonestar =  coeffvals(3);
                t1vec(i)= Tonestar*(coeffvals(2)/coeffvals(1)-1); % LL correction
                t1vec_star(i)= Tonestar;

                % Store the fitted model prediction for this voxel
                fitted_signal = abs(coeffvals(1) - coeffvals(2) * exp(-tinv / Tonestar));  % Model prediction


                % Calculate the squared residuals for each TI
                squared_residuals = (recover - fitted_signal).^2;  % Squared residuals (fit - true)^2

                % Sum the squared residuals over all inversion times and take the square root
                residuals(i) = sqrt(sum(squared_residuals));  % Residual is sqrt(sum of squared residuals)

        if t1vec(i)<0
            t1vec(i)=abs(t1vec(i));
        elseif (isnan(t1vec(i)) || isinf(t1vec(i)) || t1vec(i)>10000) % remove rubbish values, limit to 5sec max
              t1vec(i)=0; 
        end

        if t1vec_star(i)<0
            t1vec_star(i)=abs(t1vec_star(i));
        elseif (isnan(t1vec_star(i)) || isinf(t1vec_star(i)) || t1vec_star(i)>10000) % remove rubbish values, limit to 5sec max
              t1vec(i)=0; 
        end

        if isnan(residuals(i)) || isinf(residuals(i))
            residuals(i) = NaN;  % Set to 0 or NaN if invalid?
        end

    end
end

t1map = reshape(t1vec,nbrow2,nbcol2,nbslice);
t1map_star = reshape(t1vec_star,nbrow2,nbcol2,nbslice);
residuals = reshape(residuals,nbrow2,nbcol2,nbslice);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALISATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imagine(squeeze(maskedData), squeeze(registeredSeries));
% imagine(t1map, t1map_star, residuals);