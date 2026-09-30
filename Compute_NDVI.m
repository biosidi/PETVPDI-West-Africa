%% ================================================================
% MOD13A3 NDVI NORMALIZATION
% Study period: 2001-2022
%
% Purpose:
%   1. Read monthly MOD13A3 NDVI GeoTIFF files.
%   2. Apply the MOD13A3 NDVI scale factor (0.0001).
%   3. Apply the existing West Africa study-area mask.
%   4. Calculate the temporal minimum and maximum NDVI independently
%      for each pixel over the complete study period.
%   5. Normalize NDVI to the range [0,1].
%   6. Export normalized NDVI (NNDVI) as GeoTIFF files.
%
% Normalization:
%
%                 NDVI(t) - NDVImin
%   NNDVI(t) = -------------------------
%                 NDVImax - NDVImin
%
% NDVImin and NDVImax are calculated independently for each pixel
% over the complete temporal record (2001-2022).
%
% ================================================================

clear;
clc;


%% ---------------------------------------------------------------
% 1. DEFINE INPUT AND OUTPUT DIRECTORIES
% ---------------------------------------------------------------

InputDir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\2.1.West_Afr_MOD13A3_(NDVI)_2001_2022\';

OutDir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\NNDVI\';


%% ---------------------------------------------------------------
% 2. LIST ALL MOD13A3 NDVI FILES
% ---------------------------------------------------------------

DDD = dir(fullfile(InputDir, '*.tif'));

if isempty(DDD)
    error('No MOD13A3 NDVI GeoTIFF files were found.');
end

Names = {DDD.name}';
Folder = {DDD.folder}';

Filenames = cell(length(Names),1);

for i = 1:length(Names)
    Filenames{i} = fullfile(Folder{i}, Names{i});
end

fprintf('Found %d MOD13A3 NDVI files.\n', length(Filenames));


%% ---------------------------------------------------------------
% 3. LOAD THE EXISTING WEST AFRICA MASK
% ---------------------------------------------------------------
%
% RegionMask was previously generated from the West Africa
% shapefile and saved as a MAT file.
%
% Pixels inside the study region contain valid mask values,
% whereas pixels outside the study region are represented by NaN.
% ---------------------------------------------------------------

load RegionMask


%% ---------------------------------------------------------------
% 4. READ AND PREPROCESS ALL NDVI IMAGES
% ---------------------------------------------------------------
%
% MOD13A3 NDVI values are scaled using a factor of 0.0001.
%
% The West Africa mask is then applied to each monthly image.
% ---------------------------------------------------------------

IMAGES = cell(length(Filenames),1);

for k = 1:length(Filenames)

    % Read NDVI GeoTIFF
    [Image, R] = readgeoraster(Filenames{k});

    % Convert to double precision and apply scale factor
    Image = double(Image) * 0.0001;

    % Apply West Africa study-area mask
    Image = RegionMask .* Image;

    % Store processed NDVI image
    IMAGES{k,1} = Image;

end

fprintf('NDVI preprocessing completed.\n');


%% ---------------------------------------------------------------
% 5. EXTRACT DATES FROM MOD13A3 FILENAMES
% ---------------------------------------------------------------
%
% The original filenames contain the acquisition date as:
%
%       doyYYYYDDD
%
% where:
%       YYYY = year
%       DDD  = day of year
%
% Example:
%       doy2001001
%
% corresponds to 1 January 2001.
% ---------------------------------------------------------------

DateValues = zeros(length(Names),1);

for i = 1:length(Names)

    % Split original filename
    NameParts = strsplit(Names{i}, '_');

    % Preserve the original filename structure:
    % date information is contained in element 7
    DateText = NameParts{7};

    % Remove "doy"
    DateText = strrep(DateText, 'doy', '');

    % Convert YYYYDDD to numeric value
    DateValues(i) = str2double(DateText);

end


%% Convert year/day-of-year to MATLAB dates

Years = fix(DateValues / 1000);
DOY   = rem(DateValues, 1000);

Dates = datetime(Years, 1, 1) + caldays(DOY - 1);


%% ---------------------------------------------------------------
% 6. CREATE OUTPUT DIRECTORY
% ---------------------------------------------------------------

if ~exist(OutDir, 'dir')
    mkdir(OutDir);
end

fprintf('Output directory:\n%s\n', OutDir);


%% ---------------------------------------------------------------
% 7. CALCULATE PIXEL-WISE TEMPORAL MINIMUM AND MAXIMUM NDVI
% ---------------------------------------------------------------
%
% All monthly NDVI images are combined along the third dimension.
%
% For each spatial pixel:
%
%   MinNDVI = minimum NDVI during 2001-2022
%   MaxNDVI = maximum NDVI during 2001-2022
%
% Therefore, the normalization bounds are calculated independently
% for each pixel.
% ---------------------------------------------------------------

[m, n] = size(IMAGES{1});

NDVICube = reshape([IMAGES{:}], m, n, []);

% NaN-aware temporal minimum and maximum
MinNDVI = nanmin(NDVICube, [], 3);
MaxNDVI = nanmax(NDVICube, [], 3);

fprintf('Pixel-wise temporal NDVI minimum and maximum calculated.\n');


%% ---------------------------------------------------------------
% 8. CALCULATE TEMPORAL NDVI RANGE
% ---------------------------------------------------------------

NDVIRange = MaxNDVI - MinNDVI;

% Avoid division by zero where NDVImax equals NDVImin
ZeroRange = NDVIRange == 0;

NDVIRange(ZeroRange) = NaN;

fprintf('Pixels with NDVImax = NDVImin: %d\n', ...
    sum(ZeroRange(:)));


%% ---------------------------------------------------------------
% 9. NORMALIZE EACH MONTHLY NDVI IMAGE
% ---------------------------------------------------------------
%
% Pixel-wise min-max normalization:
%
%                 NDVI(t) - NDVImin
%   NNDVI(t) = -------------------------
%                 NDVImax - NDVImin
%
% The expected range of NNDVI is [0,1].
% ---------------------------------------------------------------

fprintf('\nStarting NDVI normalization...\n');

for i = 1:length(IMAGES)

    % Normalize NDVI
    NNDVI = ...
        (IMAGES{i,1} - MinNDVI) ./ NDVIRange;

    % Remove invalid values
    NNDVI(~isfinite(NNDVI)) = NaN;

    % Protect against very small floating-point rounding errors
    NNDVI(NNDVI < 0) = 0;
    NNDVI(NNDVI > 1) = 1;


    %% -----------------------------------------------------------
    % Create output date
    % -----------------------------------------------------------

    CurrentDate = Dates(i);

    YearValue  = year(CurrentDate);
    MonthValue = month(CurrentDate);


    %% -----------------------------------------------------------
    % Create output filename
    % -----------------------------------------------------------

    Outname = sprintf( ...
        'NNDVI.%04d.%02d.tif', ...
        YearValue, ...
        MonthValue);

    OutputFile = fullfile(OutDir, Outname);


    %% -----------------------------------------------------------
    % Export normalized NDVI as GeoTIFF
    % -----------------------------------------------------------

    geotiffwrite( ...
        OutputFile, ...
        NNDVI, ...
        R);


    fprintf('Normalized NDVI %d of %d: %s\n', ...
        i, length(IMAGES), Outname);

end


%% ---------------------------------------------------------------
% 10. SAVE TEMPORAL MINIMUM AND MAXIMUM NDVI
% ---------------------------------------------------------------
%
% These two GeoTIFFs document the pixel-specific normalization
% bounds used to calculate NNDVI.
% ---------------------------------------------------------------

MinFile = fullfile( ...
    OutDir, ...
    'NDVI_Temporal_Min_2001_2022.tif');

MaxFile = fullfile( ...
    OutDir, ...
    'NDVI_Temporal_Max_2001_2022.tif');


geotiffwrite( ...
    MinFile, ...
    MinNDVI, ...
    R);

geotiffwrite( ...
    MaxFile, ...
    MaxNDVI, ...
    R);


%% ---------------------------------------------------------------
% 11. FINAL CHECK
% ---------------------------------------------------------------

DDD_NNDVI = dir(fullfile(OutDir, 'NNDVI.*.tif'));

fprintf('\n============================================\n');
fprintf('NDVI NORMALIZATION COMPLETED\n');
fprintf('============================================\n');

fprintf('Normalized NNDVI files created: %d\n', ...
    length(DDD_NNDVI));

fprintf('Output folder:\n%s\n', OutDir);

fprintf('\nTemporal minimum NDVI:\n%s\n', MinFile);

fprintf('\nTemporal maximum NDVI:\n%s\n', MaxFile);

fprintf('============================================\n');