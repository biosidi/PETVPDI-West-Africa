%% GPM Precipitation Preprocessing, Resampling, and Normalization
%
% Purpose:
%   1. Read GPM precipitation data from NetCDF (.nc4) files.
%   2. Extract and mask the West African study region.
%   3. Export precipitation data as GeoTIFF files.
%   4. Resample precipitation to the MOD13A3 NDVI grid using GMT grdtrack.
%   5. Normalize precipitation independently at each pixel using the
%      temporal minimum and maximum over the study period.
%   6. Export normalized precipitation (NP) as GeoTIFF files.
%
% Normalization:
%       NP = (P - Pmin) / (Pmax - Pmin)
%
% where Pmin and Pmax are calculated independently for each pixel
% over the complete temporal record.
%
% Software:
%   MATLAB R2015
%   Generic Mapping Tools (GMT)
%
% -------------------------------------------------------------------------

clear;
clc;

%% STEP 1. Define input and output directories

% Add GMT to MATLAB path
addpath(genpath('C:\programs\gmt6\bin'));

% GPM NetCDF input directory
GPM_Input = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\GPM_ESSAI_INPUT\';

% West Africa shapefile
ShapeFile = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\WEST_AFRICA_SHAPE_ORIGINAL\WEST_AFRICA_SHAPE_ORIGINA.shp';

% NDVI input directory
NDVI_Input = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\GNDVI_INPUT\';

% Output directory for cropped GPM precipitation
GPM_TIF_Output = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\GPM_ESSAI_OUT\';

% Output directory for normalized precipitation
NP_Output = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\GPM_NP\';


%% STEP 2. List all GPM NetCDF files

DDD = dir(fullfile(GPM_Input, '*.nc4'));

Names = {DDD.name}';
Folder = {DDD.folder}';

Filenames = cell(length(Names), 1);

for i = 1:length(Names)
    Filenames{i} = fullfile(Folder{i}, Names{i});
end

% Optional:
% Display the structure and variables of the first NetCDF file.
% ncdisp(Filenames{1});


%% STEP 3. Read latitude, longitude, and precipitation

LATS = ncread(Filenames{1}, 'lat');
LONS = ncread(Filenames{1}, 'lon');

Image = ncread(Filenames{1}, 'precipitation');


%% STEP 4. Read the West Africa shapefile

[SS, ~] = shaperead(ShapeFile);


%% STEP 5. Crop the GPM grid to the West Africa bounding box

BooleanX = LONS >= min([SS.X]) & LONS <= max([SS.X]);
BooleanY = LATS >= min([SS.Y]) & LATS <= max([SS.Y]);

LONS = LONS(BooleanX);
LATS = LATS(BooleanY);


%% STEP 6. Create the West Africa spatial mask

[LonGrid, LatGrid] = meshgrid(LONS, LATS);

InMask = inpolygon( ...
    LonGrid(:), ...
    LatGrid(:), ...
    [SS.X], ...
    [SS.Y]);

Mask = reshape(InMask, size(LonGrid));

RegionMask = nan(size(Mask));
RegionMask(Mask == 1) = 1;


%% STEP 7. Apply the mask to the first precipitation image

Image = Image(BooleanY, BooleanX);
Image = RegionMask .* Image;


%% STEP 8. Create the geographic referencing object

LatLim = double([min(LATS), max(LATS)]);
LonLim = double([min(LONS), max(LONS)]);

RasterSize = size(Image);

R = georefcells( ...
    LatLim, ...
    LonLim, ...
    RasterSize, ...
    'ColumnsStartFrom', 'north');


%% STEP 9. Extract dates from GPM filenames

% The following section preserves the original filename-based
% date extraction logic.

Years  = zeros(length(Names), 1);
Months = zeros(length(Names), 1);
Days   = zeros(length(Names), 1);

Dates = cell(length(Names), 1);

for i = 1:length(Names)

    NameParts = strsplit(Names{i}, '-');
    DatePart = NameParts{2};

    DotParts = strsplit(DatePart, '.');

    DateString = DotParts{5};

    Years(i)  = str2double(DateString(1:4));
    Months(i) = str2double(DateString(5:6));
    Days(i)   = str2double(DateString(7:8));

    Dates{i} = sprintf('%04d-%02d-%02d', ...
        Years(i), Months(i), Days(i));
end


%% STEP 10. Create output directory for precipitation GeoTIFFs

if ~exist(GPM_TIF_Output, 'dir')
    mkdir(GPM_TIF_Output);
end


%% STEP 11. Crop, mask, and export each GPM image

for k = 1:length(Filenames)

    % Read precipitation
    Image = ncread(Filenames{k}, 'precipitation');

    % Crop to West Africa bounding box
    Image = Image(BooleanY, BooleanX);

    % Apply West Africa mask
    Image = RegionMask .* Image;

    % Correct raster orientation before GeoTIFF export
    Image = flipud(Image);

    % Output filename
    Outname = sprintf( ...
        'Precip.West.Africa.%02d.%04d.tif', ...
        Days(k), Months(k));

    % Write GeoTIFF
    geotiffwrite(Outname, Image, R);

    % Move output file to destination folder
    movefile(Outname, GPM_TIF_Output);

end


%% STEP 12. List the exported GPM GeoTIFF files

DDD = dir(fullfile(GPM_TIF_Output, 'Precip.West.Africa*.tif'));

Names = {DDD.name}';
Folder = {DDD.folder}';

Filenames = cell(length(Names), 1);

for i = 1:length(Names)
    Filenames{i} = fullfile(Folder{i}, Names{i});
end


%% STEP 13. Extract identifiers from precipitation filenames

FileID = cell(length(Names), 1);

for i = 1:length(Names)

    TempName = strrep(Names{i}, 'Precip.West.Africa.', '');
    TempName = strrep(TempName, '.tif', '');

    FileID{i} = TempName;

end


%% STEP 14. Read one NDVI file to define the target grid

DDD_NDVI = dir(fullfile(NDVI_Input, '*.tif'));

Names_NDVI = {DDD_NDVI.name}';
Folder_NDVI = {DDD_NDVI.folder}';

Filenames_NDVI = cell(length(Names_NDVI), 1);

for i = 1:length(Names_NDVI)
    Filenames_NDVI{i} = ...
        fullfile(Folder_NDVI{i}, Names_NDVI{i});
end

% Read NDVI grid using GMT
NDVI = gmt(['read -Tg ' Filenames_NDVI{1}]);

% Generate longitude and latitude coordinates of target NDVI grid
[NDVI_LonGrid, NDVI_LatGrid] = meshgrid(NDVI.x, NDVI.y);


%% STEP 15. Resample GPM precipitation to the NDVI grid

% The GMT grdtrack module is used to interpolate precipitation
% values at the coordinates of the NDVI grid.

IMAGES = cell(length(Filenames), 1);

[m, n] = size(NDVI_LonGrid);

for k = 1:length(Filenames)

    GMT_Command = ['grdtrack -N -G' Filenames{k}];

    Grid_Positions = gmt( ...
        GMT_Command, ...
        [NDVI_LonGrid(:), NDVI_LatGrid(:)]);

    % Third column contains interpolated precipitation values
    Grid_Positions = Grid_Positions.data(:,3);

    % Convert interpolated values back to raster format
    Resampled_Precip = reshape(Grid_Positions, m, n);

    % Store each resampled precipitation raster
    IMAGES{k,1} = Resampled_Precip;

end


%% STEP 16. Create output directory for normalized precipitation

if ~exist(NP_Output, 'dir')
    mkdir(NP_Output);
end


%% STEP 17. Create geographic reference for the NDVI target grid

LatLim = [min(NDVI.y), max(NDVI.y)];
LonLim = [min(NDVI.x), max(NDVI.x)];

RasterSize = size(NDVI.z);

R = georefcells( ...
    LatLim, ...
    LonLim, ...
    RasterSize, ...
    'ColumnsStartFrom', 'north');


%% STEP 18. Calculate pixel-wise temporal minimum and maximum precipitation

[m, n] = size(NDVI.z);

% Convert the cell array into a 3-D precipitation data cube:
% rows x columns x time

PrecipCube = reshape([IMAGES{:}], m, n, []);

% MATLAB R2015-compatible NaN-aware temporal minimum and maximum.
% The Statistics and Machine Learning Toolbox provides nanmin/nanmax.

MinP = nanmin(PrecipCube, [], 3);
MaxP = nanmax(PrecipCube, [], 3);


%% STEP 19. Normalize precipitation

% Pixel-wise min-max normalization:
%
%               P(t) - Pmin
%       NP(t) = -------------
%               Pmax - Pmin
%
% Pmin and Pmax are calculated independently for each pixel
% across the complete temporal record.

PPP = parpool(8);

parfor i = 1:length(IMAGES)

    NPrecip = ...
        (IMAGES{i,1} - MinP) ./ ...
        (MaxP - MinP);

    % Output filename
    Outname = ['NP_GPM.' FileID{i} '.tif'];

    % Export normalized precipitation
    geotiffwrite(Outname, flipud(NPrecip), R);

    % Move output to destination folder
    movefile(Outname, NP_Output);

end

delete(PPP);


%% END OF SCRIPT