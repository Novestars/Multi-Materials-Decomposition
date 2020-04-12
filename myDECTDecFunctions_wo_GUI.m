InputImageE1_path = 'D:\Backup\Desktop\vessel12\DECTDec-master\Test Data\Test Data\Mouse head specimen\80kVp 16bit unsigned\Mousehead_80kV_FOV2_HU_0620.tif';
InputImageE2_path = strrep(InputImageE1_path, '80kV', '40kV');

InputImageE1 = imread(InputImageE1_path);
InputImageE2 = imread(InputImageE2_path);


setting_E1 = load('D:\Backup\Desktop\vessel12\DECTDec-master\Test Data\Test Data\Mouse head specimen\Mouse_head_80kVp.mat');
setting_E2 = load('D:\Backup\Desktop\vessel12\DECTDec-master\Test Data\Test Data\Mouse head specimen\Mouse_head_40kVp.mat');

E1 = [setting_E1.EnergySetting.Mat1(1), setting_E1.EnergySetting.Mat2(1),setting_E1.EnergySetting.Mat3(1)];
E2 = [setting_E2.EnergySetting.Mat1(1), setting_E2.EnergySetting.Mat2(1),setting_E2.EnergySetting.Mat3(1)];

[OutputImage1,OutputImage2,OutputImage3] = ProcessImages(InputImageE1, InputImageE2, E1, E2);

figure
subplot(2,3,1)
imagesc(InputImageE1);
subplot(2,3,2)
imagesc(InputImageE2);
subplot(2,3,4)
imagesc(OutputImage1);
subplot(2,3,5)
imagesc(OutputImage2);
subplot(2,3,6)
imagesc(OutputImage3);

function [OutputImage1,OutputImage2,OutputImage3] = ProcessImages(InputImageE1,InputImageE2,Energy1, Energy2)

    sizeInE1 = size(InputImageE1(:,:,1));
    sizeInE2 = size(InputImageE2(:,:,1));
    if sizeInE1 ~= sizeInE2
        return 
    end

    % Preparing the array with absorption properties
    E1Mat1 = Energy1(1);            % Mat1 = absorption of material 1 at energy 1
    E1Mat2 = Energy1(2);            % Mat2 = absorption of material 2 at energy 1
    E1Mat3 = Energy1(3);            % Mat3 = absorption of material 3 at energy 1
    E2Mat1 = Energy2(1);            % Mat1 = absorption of material 1 at energy 2
    E2Mat2 = Energy2(2);            % Mat2 = absorption of material 2 at energy 2
    E2Mat3 = Energy2(3);            % Mat3 = absorption of material 3 at energy 2

    A = [E1Mat1 E1Mat2 E1Mat3; E2Mat1 E2Mat2 E2Mat3; 1 1 1];
    AINV = inv(A);

    % Reshape Image 1 to single row vector
    vector1 = reshape(InputImageE1, [], 1);
    vector1rot = rot90(vector1);

    % Reshape Image 2 to single row vector
    vector2 = reshape(InputImageE2, [], 1);
    vector2rot = rot90(vector2);

    % Assemble matrix (3 rows; Enery 1, Energy 2, 1)
    matrix = vertcat(vector1rot,vector2rot);
    matrix(3,:) = 1;
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

    %Extract material 3 and reshape matrix to image
    Material3 = matrixbackrot(:,3);
    decomposite3 = reshape(Material3, [sizeInE1, 1]);
    decomposite3 = decomposite3 * 65535;
    OutputImage3 = uint16(decomposite3);

end % function ProcessImage





