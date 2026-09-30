%% ================================================================
% PETVPDI CALCULATION
% Study period: 2001-2022
%
% Purpose:
%   1. Read normalized PET (NPET), NDVI (NNDVI), and
%      precipitation (NP) GeoTIFF files.
%   2. Verify that the three datasets contain the same number
%      of monthly observations.
%   3. Calculate PETVPDI for each month.
%   4. Export monthly PETVPDI as GeoTIFF files.
%
% PETVPDI formulation:
%
%   PETVPDI = sqrt[(1 - NPET)^2 + NNDVI^2 + NP^2]
%
% The index represents the Euclidean distance from the normalized
% dry reference state D = (1,0,0).
%
% The theoretical PETVPDI range is:
%
%       0 <= PETVPDI <= sqrt(3)
%
% Lower values indicate conditions closer to the normalized dry
% reference state, whereas higher values indicate conditions
% farther from the dry reference state.
%
% ================================================================

clear;
clc;


%% ---------------------------------------------------------------
% 1. DEFINE INPUT DIRECTORIES
% ---------------------------------------------------------------

NPET_Dir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\NPET\NPET\';

NNDVI_Dir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\NNDVI\NNDVI\';

NP_Dir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\NP_GPM\NP\';


%% ---------------------------------------------------------------
% 2. DEFINE OUTPUT DIRECTORY
% ---------------------------------------------------------------

OutDir = ...
    'G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\PETVPDI_GPM\';

if ~exist(OutDir, 'dir')
    mkdir(OutDir);
end


%% ---------------------------------------------------------------
% 3. LIST NORMALIZED PET FILES
% ---------------------------------------------------------------

DDD_NPET = dir(fullfile(NPET_Dir, '*.tif'));

if isempty(DDD_NPET)
    error('No normalized PET files were found.');
end

Names_NPET = {DDD_NPET.name}';
Folder_NPET = {DDD_NPET.folder}';

Filenames_NPET = cell(length(Names_NPET),1);

for i = 1:length(Names_NPET)
    Filenames_NPET{i} = ...
        fullfile(Folder_NPET{i}, Names_NPET{i});
end


%% ---------------------------------------------------------------
% 4. LIST NORMALIZED NDVI FILES
% ---------------------------------------------------------------

DDD_NNDVI = dir(fullfile(NNDVI_Dir, '*.tif'));

if isempty(DDD_NNDVI)
    error('No normalized NDVI files were found.');
end

Names_NNDVI = {DDD_NNDVI.name}';
Folder_NNDVI = {DDD_NNDVI.folder}';

Filenames_NNDVI = cell(length(Names_NNDVI),1);

for i = 1:length(Names_NNDVI)
    Filenames_NNDVI{i} = ...
        fullfile(Folder_NNDVI{i}, Names_NNDVI{i});
end


%% ---------------------------------------------------------------
% 5. LIST NORMALIZED PRECIPITATION FILES
% ---------------------------------------------------------------

DDD_NP = dir(fullfile(NP_Dir, '*.tif'));

if isempty(DDD_NP)
    error('No normalized precipitation files were found.');
end

Names_NP = {DDD_NP.name}';
Folder_NP = {DDD_NP.folder}';

Filenames_NP = cell(length(Names_NP),1);

for i = 1:length(Names_NP)
    Filenames_NP{i} = ...
        fullfile(Folder_NP{i}, Names_NP{i});
end


%% ---------------------------------------------------------------
% 6. VERIFY NUMBER OF INPUT FILES
% ---------------------------------------------------------------

nNPET  = length(Filenames_NPET);
nNNDVI = length(Filenames_NNDVI);
nNP    = length(Filenames_NP);

fprintf('NPET files  : %d\n', nNPET);
fprintf('NNDVI files : %d\n', nNNDVI);
fprintf('NP files    : %d\n', nNP);

assert( ...
    nNPET == nNNDVI && nNPET == nNP, ...
    'NPET, NNDVI, and NP must contain the same number of files.');

fprintf('All three datasets contain the same number of files.\n');


%% ---------------------------------------------------------------
% 7. EXTRACT DATE INFORMATION
% ---------------------------------------------------------------
%
% This section follows the naming convention of the normalized
% precipitation files:
%
%       NP.YYYY.MM.tif
%
% If your actual filenames use another convention, this section
% should be adjusted accordingly.
% ---------------------------------------------------------------

DateX = cell(nNP,1);

for i = 1:nNP

    CurrentName = Names_NP{i};

    % Remove prefix
    CurrentName = strrep(CurrentName, 'NP.', '');

    % Remove file extension
    CurrentName = strrep(CurrentName, '.tif', '');

    DateX{i} = CurrentName;

end


%% ---------------------------------------------------------------
% 8. READ GEOGRAPHIC REFERENCE FROM FIRST NPET FILE
% ---------------------------------------------------------------

[FirstNPET, R] = readgeoraster(Filenames_NPET{1});

FirstNPET = double(FirstNPET);

[nRows,nCols] = size(FirstNPET);

fprintf('Raster dimensions: %d x %d\n', ...
    nRows,nCols);


%% ---------------------------------------------------------------
% 9. CALCULATE PETVPDI FOR EACH MONTH
% ---------------------------------------------------------------
%
% The three input variables have already been normalized:
%
%       NPET  = normalized PET
%       NNDVI = normalized NDVI
%       NP    = normalized precipitation
%
% Therefore:
%
% PETVPDI = sqrt[(1-NPET)^2 + NNDVI^2 + NP^2]
%
% No additional normalization is performed in this step.
% ---------------------------------------------------------------

fprintf('\n============================================\n');
fprintf('STARTING PETVPDI CALCULATION\n');
fprintf('============================================\n');


for i = 1:nNPET

    %% Read normalized PET
    [NPET, R_NPET] = ...
        readgeoraster(Filenames_NPET{i});

    NPET = double(NPET);


    %% Read normalized NDVI
    [NNDVI, R_NNDVI] = ...
        readgeoraster(Filenames_NNDVI{i});

    NNDVI = double(NNDVI);


    %% Read normalized precipitation
    [NP, R_NP] = ...
        readgeoraster(Filenames_NP{i});

    NP = double(NP);


    %% -----------------------------------------------------------
    % Verify raster dimensions
    % -----------------------------------------------------------

    if ~isequal(size(NPET),size(NNDVI),size(NP))

        error( ...
            'Raster dimensions do not match for observation %d.', ...
            i);

    end


    %% -----------------------------------------------------------
    % Remove non-finite values
    % -----------------------------------------------------------

    NPET(~isfinite(NPET))   = NaN;
    NNDVI(~isfinite(NNDVI)) = NaN;
    NP(~isfinite(NP))       = NaN;


    %% -----------------------------------------------------------
    % Calculate PETVPDI
    % -----------------------------------------------------------

    PETVPDI = sqrt( ...
        (1 - NPET).^2 + ...
        NNDVI.^2 + ...
        NP.^2 );


    %% -----------------------------------------------------------
    % Remove invalid output
    % -----------------------------------------------------------

    PETVPDI(~isfinite(PETVPDI)) = NaN;


    %% -----------------------------------------------------------
    % Create output filename
    % -----------------------------------------------------------

    Outname = ...
        ['PETVPDI_GPM.' DateX{i} '.tif'];

    OutputFile = ...
        fullfile(OutDir, Outname);


    %% -----------------------------------------------------------
    % Export PETVPDI as GeoTIFF
    % -----------------------------------------------------------

    geotiffwrite( ...
        OutputFile, ...
        PETVPDI, ...
        R_NPET);


    fprintf( ...
        'PETVPDI %d of %d: %s\n', ...
        i, nNPET, Outname);


    %% Free temporary variables
    clear NPET NNDVI NP PETVPDI

end


%% ---------------------------------------------------------------
% 10. FINAL CHECK
% ---------------------------------------------------------------

DDD_Output = ...
    dir(fullfile(OutDir, 'PETVPDI_GPM.*.tif'));

fprintf('\n============================================\n');
fprintf('PETVPDI CALCULATION COMPLETED\n');
fprintf('============================================\n');

fprintf('PETVPDI files created: %d\n', ...
    length(DDD_Output));

fprintf('Output folder:\n%s\n', ...
    OutDir);

fprintf('============================================\n');