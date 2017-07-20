function [SFNR, SNR] = roi_snr(test_image, spot_focus, spot_r, nDisdaqs)
%function [SFNR, SNR] = roi_snr(test_image, spot_focus, spot_r, [nDisdaqs])
% test_image = full path to the 4D image of interest 
% spot_focus = 3D voxel space input to center sphere on [x, y, z] 1X3
% spot_r = radius of spotlight sphere in mm
%    (only voxels within the image indices are used)
% nDisdaqs = number of volumes to remove from start of run before calculating SNR/SFNR
%%% Version 1.0  5/10/08 - JJM
%%% Version 1.1  6/08/08 - RMC: functioned, Modified to make and use a sphere ROI
%%% Version 1.2  7/24/08 - JJM: option to remove DisDaqs, removes dead voxels from ROIs
%%% Reads in an image file, creates masks of specified ROI's and generates
%%% SNR and SFNR values for each ROI




a = readmr(test_image);       % loads the MR image we want. 

if nargin == 3
    a_img_dims = [a.info.dimensions.size];
    a_vox_dims = [a.info.dimensions.spacing];
    n_vox_r = floor(spot_r ./ a_vox_dims(1:3));
    timepoints = 1:length(a.data(1,1,1,:));
    n_timepoints = length(timepoints);

elseif nargin > 3                                   % block of code to remove disdaqs, if specified
    a.data = a.data(:,:,:,(nDisdaqs + 1):end);
    a.info.dimensions(4).size = (a.info.dimensions(4).size - nDisdaqs);
    
    a_img_dims = [a.info.dimensions.size];           
    a_vox_dims = [a.info.dimensions.spacing];
    n_vox_r = floor(spot_r ./ a_vox_dims(1:3));
    timepoints = 1:length(a.data(1,1,1,:));
    n_timepoints = length(timepoints);
    
end

%make mask starter
mask = zeros(size(a.data));             % creates an empty mask the size of the bxh image we're dealing with
mask = mask(:,:,:,1);                   % collapse into 3-dimensions


%create sphere within given radius and within image bounds
sphere_voxels = [];
for z = (-1*n_vox_r(3)):1:n_vox_r(3)    %for each z value
    if (z+spot_focus(3)>0) & (z+spot_focus(3)<=a_img_dims(3))  %if within z image dim
        for y = (-1*n_vox_r(2)):1:n_vox_r(2)    %for each y value
            if (y+spot_focus(2)>0) & (y+spot_focus(2)<=a_img_dims(2))   %if within y image dim
                for x = (-1*n_vox_r(1)):1:n_vox_r(1)  %for each x value
                    if (x+spot_focus(1)>0) & (x+spot_focus(1)<=a_img_dims(1))  %two qualifications %1. inside the image, checked y and z above
                        if pdist([spot_focus; [x*a_vox_dims(1) y*a_vox_dims(2) z*a_vox_dims(3)]+spot_focus], 'euclidean') <= spot_r  %2. inside spot_r
                            mask(spot_focus(1)+x, spot_focus(2)+y, spot_focus(3)+z) = 100;    %set mask voxel value to something other than 0
                            %sphere_voxels = [sphere_voxels; spot_focus + [x y z]]    %write the sphere voxel to the valid list
                        end   %end if within radius
                    end   %end if inside x image dims
                end  %end for each x
            end   %end if inside y image dims
         end   %end for each y
    end  % end if inside z img dims
end   %end for each z


%tester
%mask(:,:,:) = 0;                        % reset mask
%mask(38:48, 39:45, 1:8) = 100;


%%% Extract ROI's from masks

I = find(mask>0);                   % 'I' will be a linear index of the non-zero entities in mask
[x,y,z] = ind2sub(size(a.data), I);     % Converts the linear index into the 3-d coordinates of each point

roi = [];                           
for j = 1:length(a.data(1,1,1,:))       % 'j' will be 1:number of time points
    tmp = [];                           % temporary holding matrix
    for i = 1:length(x)                 % 'i' will be the number of voxels in the masked section
        tmp = [tmp; a.data(x(i),y(i),z(i),j)];  % tmp will be the intensity of each voxel in the masked section at timepoint 'j'
    end
     roi = [roi tmp];            % roi will be a 2D matrix of voxel intensities (rows) across each time point (columns)        
end

%%% Remove dead voxels from masked region (when voxel intensity = 0 for all time points)

vox_intensity_sum = sum(roi,2);     % Sum across time for each voxel's intensity
dead_vox = find(vox_intensity_sum == 0);  % row # of dead voxels
for c = 1:length(dead_vox)                  
    roi(dead_vox(c),:) = [];              % deletes the rows in "roi" where the dead voxels are
end
c = [];                             % reset c




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Calculate SNR and SFNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SNR and SFNR are being calculated using methods outlined in Friedman &
% Glover (2006), Report on a Multicenter fMRI Quality Assurance Protocol,
% Journal of Magnetic Resonance Imaging, 23:827-839

%%%%%%%% Signal-to-Fluctuation-Noise Ratio (SFNR) %%%%%%%%%%

%%%% Signal Image: The average of each voxel in the sample across all time points

roi_signal_image = mean(roi, 2);     % array of average vox intensity across time for each vox in the masked region



%%%% Temporal Fluctuation Noise Image: Std. deviation of each voxel in the time series after detrending with 2nd-order polynomial
detrended_roi = [];                % empty matrix to write all of the detrended voxels into

for c = 1:length(roi(:,1))
    voxel = roi(c,:);            % will go through each voxel in the sample and detrend its intensity over time
    quad_fit = polyfit(timepoints, voxel, 2);       % fits a 2nd-order polynomial to the data
    trend = polyval(quad_fit, timepoints);          
    resid = voxel - trend;                          % subtracts the trend line from the orig data
    detrended_roi = [detrended_roi; resid];
end

roi_tempfluc_image = std(detrended_roi,0,2);        % will calculate the std deviation for each voxel across time



%%%% Generate the SFNR Image and Summary value
roi_sfnr_image = roi_signal_image./roi_tempfluc_image;  % array with individual sfnr for each voxel
SFNR = mean(roi_sfnr_image);                             % average sfnr across all voxels 


%%%%%%%% Signal-to-Noise Ratio (SNR) %%%%%%%%%%
% roi sphere mask
 roi_odd = sum(roi(:,1:2:end),2);   % Sum of odd images (sum of each voxel's intensity over time)
 roi_even = sum(roi(:,2:2:end),2);  % Sum of even images
 roi_diff = roi_odd - roi_even;       % difference image - raw measure of spatial noise
 
 roi_variance_summary = var(roi_diff); % the amount of variance in the spatial noise difference image
 roi_sum_image = sum(roi,2);   % sum of all images (sum of each voxel's intensity over time)
 roi_avg_image = roi_sum_image/n_timepoints; 
 
 SNR = (mean(roi_avg_image)/sqrt(roi_variance_summary/n_timepoints));

 
 

