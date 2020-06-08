% This function is used to investigate the performance of the dark blood
% algorithm on the aorty. The two metrial decomposition method used here is tranferred from the three material
% decomposition method. 

%load label

% Load Data 
dir_path = 'D:\process_XZ\MLIII_Contrast_reconstructions';
subject_name = 'H-1489_MESAL-8012091_Mono';
suf_40 = '40.nii';
suf_70 = '70.nii';

InputImageE1_path = fullfile(dir_path,strcat( subject_name, suf_40));
InputImageE2_path = fullfile(dir_path,strcat( subject_name, suf_70));

InputImageE13D = permute( niftiread(InputImageE1_path),[3,2,1]);
InputImageE23D = permute( niftiread(InputImageE2_path),[3,2,1]);


%processing

% mask_vessel_E1 = hh.*double(InputImageE13D);
% mask_vessel_E2 = hh.*double(InputImageE23D);
% mask_vessel_E1(mask_vessel_E1==0)=-1000;
% mask_vessel_E2(mask_vessel_E2==0)=-1000;

% mean_ct_seg_E1 =  zeros(size(segs2,1),1);
% mean_ct_seg_E2 =  zeros(size(segs2,1),1);
% 
% 
% for k = 1:length(segs2)
%     unrotind = round(segs2{k}+1);
%     unrotind = unrotind(:,[2,1,3]);
%     my_ind = sub2ind(size(roi_value),unrotind(:,1),unrotind(:,2),unrotind(:,3));
%     mean_ct_seg_E1(k) = median(InputImageE13D(my_ind)); % One is added since coordinates are zero based
%     
% 
%     mean_ct_seg_E2(k) = median(InputImageE23D(my_ind)); % One is added since coordinates are zero based
%    
%  end


% Filtering
%InputImageE13D = imgaussfilt3(InputImageE13D,1);
%InputImageE23D = imgaussfilt3(InputImageE23D,1);
% a = reshape(arr,512,3008)';
% wn = reshape(E_Temp2',512,47,64);

% E_Temp2 	=   WNNM( a, 1, 3, ones(512,3008)', 100); % WNNM Estimation

% Prepocess Data: slice selection and conversion bewteen HU and linear
% attenuation
water_mu_80kev = 7.865E-02;
water_mu_150kev = 5.754E-02;
water_mu_40kev = 1.061E-01;
water_mu_70kev = (8.956E-02+7.865E-02)/2;
air_mu_80kev = 7.074E-02;
air_mu_150kev = 5.175E-02;
air_mu_40kev = 9.549E-02;
air_mu_70kev = (8.055E-02+7.074E-02)/2;


shape = size(InputImageE13D);

dc1 = zeros(shape)-1000;
dc2 = zeros(shape)-1000;
dc3 = zeros(shape)-1000;

mask_vessel_E1 = HU2LinearAtten(InputImageE13D, water_mu_40kev, air_mu_40kev);
mask_vessel_E2 = HU2LinearAtten(InputImageE23D, water_mu_70kev, air_mu_70kev);

density_water = 0.99802;
density_air = 1.2041E-3;
density_fat = 0.9094 ;


for i = 1:size(segs2,1)
    

    % % Energy setting

%     ROI_E1 = mean_ct_seg_E1(i);
%     ROI_E2 = mean_ct_seg_E2(i);

%     mu_ROI_40kev = (500/1000)*(water_mu_40kev*density_water - air_mu_40kev*density_air)+water_mu_40kev*density_water;
%     mu_ROI_70kev = (173/1000)*(water_mu_70kev*density_water - air_mu_70kev*density_air)+water_mu_70kev*density_water;

%     mu_ROI2_40kev = (-14/1000)*(water_mu_40kev*density_water - air_mu_40kev*density_air)+water_mu_40kev*density_water;
%     mu_ROI2_70kev = (-28/1000)*(water_mu_70kev*density_water - air_mu_70kev*density_air)+water_mu_70kev*density_water;
    mask = find(ct_seg ==i);
    if isempty(mask) 
        continue;
    end
    
    InputVesselE1 = mask_vessel_E1(mask);
    InputVesselE2 = mask_vessel_E2(mask);
    
    
    E1 = [median(InputVesselE1), water_mu_40kev*density_water];
    E2 = [median(InputVesselE2), water_mu_70kev*density_water];

    [OutputImage1,OutputImage2] = ProcessImages(InputVesselE1, InputVesselE2, E1, E2);
    dc1(mask) = OutputImage1;
    dc2(mask) = OutputImage2;

end

PlotImages(InputImageE13D, InputImageE23D, dc1,dc2,hh, 160);
dc1(dc1>10) = 10;
dc1(dc1<-10) = -10;
dc2(dc1>10) = 10;
dc2(dc1<-10) = -10;

 niftiwrite(int16(rescale( dc1)*65535),'result/dc1.nii')
 niftiwrite(int16(rescale( dc2)*65535),'result/dc2.nii')


function  PlotImages(InputImageE13D, InputImageE23D, InputImage1,InputImage2,InputImage3, slice)

figure
subplot(2,3,1)
imagesc(InputImageE13D(:,:,slice));
colormap gray
subplot(2,3,2)
imagesc(InputImageE23D(:,:,slice));
subplot(2,3,4)
imagesc(InputImage1(:,:,slice));
subplot(2,3,5)
imagesc(InputImage2(:,:,slice));
subplot(2,3,6)
imagesc(InputImage3(:,:,slice));
end

function [OutputImage1,OutputImage2] = ProcessImages(InputImageE1,InputImageE2,Energy1, Energy2)

    sizeInE1 = size(InputImageE1(:,:,1));
    sizeInE2 = size(InputImageE2(:,:,1));
    if sizeInE1 ~= sizeInE2
        return 
    end

    % Preparing the array with absorption properties
    E1Mat1 = Energy1(1);            % Mat1 = absorption of material 1 at energy 1
    E1Mat2 = Energy1(2);            % Mat2 = absorption of material 2 at energy 1
    E2Mat1 = Energy2(1);            % Mat1 = absorption of material 1 at energy 2
    E2Mat2 = Energy2(2);            % Mat2 = absorption of material 2 at energy 2

    A = [E1Mat1 E1Mat2 ; E2Mat1 E2Mat2];
    AINV = inv(A);

    % Reshape Image 1 to single row vector
    vector1 = reshape(InputImageE1, [], 1);
    vector1rot = rot90(vector1);

    % Reshape Image 2 to single row vector
    vector2 = reshape(InputImageE2, [], 1);
    vector2rot = rot90(vector2);

    % Assemble matrix (3 rows; Enery 1, Energy 2, 1)
    matrix = vertcat(vector1rot,vector2rot);
    matrixdouble = double(matrix);

    % Dual energy decomposition
    Result(:,:) = AINV * matrixdouble(:,:);

    % Build matrix from images
    Resultflip = flipud(Result);
    matrixbackrot = rot90(Resultflip, -1);

    %Extract material 1 and reshape matrix to image
    Material1 = matrixbackrot(:,1);
    decomposite1 = reshape(Material1, [sizeInE1, 1]);
    %decomposite1 = decomposite1 * 65535;
    OutputImage1 = double(decomposite1);

    %Extract material 2 and reshape matrix to image
    Material2 = matrixbackrot(:,2);
    decomposite2 = reshape(Material2, [sizeInE1, 1]);
    %decomposite2 = decomposite2 * 65535;
    OutputImage2 = double(decomposite2);

end % function ProcessImage


function Output=HU2LinearAtten(Input, mass_attenuation_water, mass_attenuation_air)

density_water = 0.99802;
density_air = 1.2041E-3;
linear_attenuation_water = mass_attenuation_water * density_water;
linear_attenuation_air = mass_attenuation_air * density_air;

Output = (double(Input)./1000).*(linear_attenuation_water-linear_attenuation_air) + linear_attenuation_water;
end


