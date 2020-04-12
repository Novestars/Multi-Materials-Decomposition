% This function is used to investigate the performance of the dark blood
% algorithm on the aorty. The two metrial decomposition method used here is tranferred from the three material
% decomposition method. The two candidate material are ROI/water. 

% Load Data 
dir_path = 'D:\process_XZ\MLIII_Contrast_raw_CT';
subject_name = 'H-1489_MESAL-8012091_MESAL-8012091-E3_FRC-DE-0.75-Qr40d-5-DET';
suf_80 = '-A-80kV_54882745.nii';
suf_150 = '-B-Sn150kV_54883253.nii';

InputImageE1_path = fullfile(dir_path,strcat( subject_name, suf_80));
InputImageE2_path = fullfile(dir_path,strcat( subject_name, suf_150));

InputImageE13D = niftiread(InputImageE1_path);
InputImageE23D = niftiread(InputImageE2_path);

% Filtering
InputImageE13D = imgaussfilt3(InputImageE13D,1);
InputImageE23D = imgaussfilt3(InputImageE23D,1);
% 
% Prepocess Data: slice selection and conversion bewteen HU and linear
% attenuation
water_mu_80kev = 7.865E-02;
water_mu_150kev = 5.754E-02;
air_mu_80kev = 7.074E-02;
air_mu_150kev = 5.175E-02;

shape = size(InputImageE13D);

dc1 = zeros(shape);
dc2 = zeros(shape);
dc3 = zeros(shape);

for i = 160%1:shape(3)
    InputImageE1=InputImageE13D(:,:,i);
    InputImageE2=InputImageE23D(:,:,i);

    InputImageE1 = HU2LinearAtten(InputImageE1, water_mu_80kev, air_mu_80kev);
    InputImageE2 = HU2LinearAtten(InputImageE2, water_mu_150kev, air_mu_150kev);

    % % Energy setting
    density_water = 0.99802;
    density_air = 1.2041E-3;
    density_fat = 0.9094 ;

    mu_ROI_80kev = (229/1000)*(water_mu_80kev*density_water - air_mu_80kev*density_air)+water_mu_80kev*density_water;
    mu_ROI_150kev = (93/1000)*(water_mu_150kev*density_water - air_mu_150kev*density_air)+water_mu_150kev*density_water;

    mu_ROI2_80kev = (-14/1000)*(water_mu_80kev*density_water - air_mu_80kev*density_air)+water_mu_80kev*density_water;
    mu_ROI2_150kev = (-28/1000)*(water_mu_150kev*density_water - air_mu_150kev*density_air)+water_mu_150kev*density_water;


    E1 = [mu_ROI_80kev, water_mu_80kev*density_water];
    E2 = [mu_ROI_150kev, water_mu_150kev*density_water];

    [OutputImage1,OutputImage2] = ProcessImages(InputImageE1, InputImageE2, E1, E2);
    dc1(:,:,i) = OutputImage1;
    dc2(:,:,i) = OutputImage2;

end
PlotImages(InputImageE13D, InputImageE23D, dc1,dc2,dc3, 160);

% niftiwrite(int16(dc1),'result/dc1.nii')
% niftiwrite(int16(dc2),'result/dc2.nii')
% niftiwrite(int16(dc3),'result/dc3.nii')

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
    decomposite1 = decomposite1 * 65535;
    OutputImage1 = uint16(decomposite1);

    %Extract material 2 and reshape matrix to image
    Material2 = matrixbackrot(:,2);
    decomposite2 = reshape(Material2, [sizeInE1, 1]);
    decomposite2 = decomposite2 * 65535;
    OutputImage2 = uint16(decomposite2);

end % function ProcessImage


function Output=HU2LinearAtten(Input, mass_attenuation_water, mass_attenuation_air)

density_water = 0.99802;
density_air = 1.2041E-3;
linear_attenuation_water = mass_attenuation_water * density_water;
linear_attenuation_air = mass_attenuation_air * density_air;

Output = (double(Input)./1000).*(linear_attenuation_water-linear_attenuation_air) + linear_attenuation_water;
end


