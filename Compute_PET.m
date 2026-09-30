%% ================================================================
% MOD16A2GF PET: 8-DAY -> MONTHLY -> Resample to 1-km ANALYSIS GRID
% Study period: 2001-2022
%
% IMPORTANT:
% MOD16A2GF PET:
%   - Original unit: kg m^-2 per composite
%   - PET represents an 8-day TOTAL
%   - Scale factor: 0.1
%   - 1 kg m^-2 water = 1 mm water
%   - Special/fill codes: 32761-32767
%
% Processing:
%   1. Read one 8-day raster at a time
%   2. Remove fill/special values BEFORE scaling
%   3. Apply scale factor 0.1
%   4. Allocate each composite to calendar months according
%      to temporal overlap
%   5. Sum contributions -> monthly PET (mm month^-1)
%   6. Mask outside West Africa as NaN
%   7. Resample to the target ~1-km grid
%
% Memory-efficient: no 264-month 3-D data cube.
% ================================================================

clear;
clc;

% 1. PATHS
% Add GMT to MATLAB path
addpath(genpath("C:\programs\gmt6\bin"))
InputDir = ...
"G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\7.1.West_Afr_MOD16A2GF_(PET)_2001_2022";

TemporaryOutDir = ...
"G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\PET.Monthly.Corrected";

%%% Open one of the NDVI files:
DDD_NDVI = dir("G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\West_Afr_MOD13A3_NDVI_2001_2022\*.tif");
Names_NDVI = {DDD_NDVI.name}';
Names_NDVI = string(Names_NDVI);
Folder_NDVI = {DDD_NDVI.folder}';
Folder_NDVI = string(Folder_NDVI);
Filenames_NDVI = fullfile(Folder_NDVI, Names_NDVI);

NDVI = gmt(char("read -Tg " + Filenames_NDVI(1)));
[NDVI_LonGrid, NDVI_LatGrid] = meshgrid(NDVI.x, NDVI.y);

LatLim = [min(NDVI.y), max(NDVI.y)];
LonLim = [min(NDVI.x), max(NDVI.x)];
RasterSize = size(NDVI.z);
R = georefcells(LatLim, LonLim, RasterSize, "ColumnsStartFrom","north");


OutDir = "G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\PET.Monthly.Resampled.Corrected";


if ~exist(TemporaryOutDir,'dir')
    mkdir(TemporaryOutDir);
end

if ~exist(OutDir,'dir')
    mkdir(OutDir);
end


%% ---------------------------------------------------------------
% 2. LOAD EXISTING WEST AFRICA MASK / TARGET GRID INFORMATION
% ---------------------------------------------------------------

load RegionMask_PET

% This assumes your existing workspace/MAT file also provides:
%
% NDVI_LonGrid
% NDVI_LatGrid
% m
% n
%
% exactly as in your original code.
%
% If these variables are stored elsewhere, load that MAT file here.


%% ---------------------------------------------------------------
% 3. LIST ORIGINAL 8-DAY MOD16A2GF PET FILES
% ---------------------------------------------------------------

DDD_PET = dir(fullfile(InputDir,'*.tif'));

if isempty(DDD_PET)
    error('No PET GeoTIFF files found.');
end

Names_PET = string({DDD_PET.name})';
Folder_PET = string({DDD_PET.folder})';
Filenames_PET = fullfile(Folder_PET,Names_PET);

fprintf('Found %d original PET files.\n',numel(Filenames_PET));


%% ---------------------------------------------------------------
% 4. EXTRACT START DATE FROM "doyYYYYDDD"
% ---------------------------------------------------------------

startDate = NaT(numel(Names_PET),1);

for i = 1:numel(Names_PET)

    token = regexp(Names_PET(i), ...
        'doy(\d{4})(\d{3})', ...
        'tokens','once');

    if isempty(token)
        error('Cannot extract date from filename: %s',Names_PET(i));
    end

    yr  = str2double(token{1});
    doy = str2double(token{2});

    startDate(i) = datetime(yr,1,1) + days(doy-1);

end


%% Sort chronologically

[startDate,sortIndex] = sort(startDate);

Names_PET = Names_PET(sortIndex);
Filenames_PET = Filenames_PET(sortIndex);


%% Restrict study period

keep = year(startDate) >= 2001 & year(startDate) <= 2022;

startDate = startDate(keep);
Names_PET = Names_PET(keep);
Filenames_PET = Filenames_PET(keep);


%% ---------------------------------------------------------------
% 5. DETERMINE COMPOSITE END DATES
% ---------------------------------------------------------------
%
% MOD16 composites normally represent 8-day totals.
%
% Using the next observation date is safer than blindly assuming
% eight days for the final composite of each year.
% ---------------------------------------------------------------

nInput = numel(startDate);

endDate = NaT(nInput,1);

for i = 1:nInput

    yr = year(startDate(i));

    if i < nInput && year(startDate(i+1)) == yr

        % Day before next composite
        endDate(i) = startDate(i+1) - days(1);

    else

        % Final composite ends on December 31
        endDate(i) = datetime(yr,12,31);

    end

end

nDaysComposite = days(endDate-startDate) + 1;


%% Check durations

fprintf('\nComposite duration summary:\n');
disp(unique(nDaysComposite));


%% ---------------------------------------------------------------
% 6. DEFINE 264 CALENDAR MONTHS
% ---------------------------------------------------------------

monthDates = ...
    (datetime(2001,1,1):calmonths(1):datetime(2022,12,1))';

nMonths = numel(monthDates);

assert(nMonths == 264,'Expected 264 months.');


%% ---------------------------------------------------------------
% 7. READ FIRST FILE FOR DIMENSIONS / REFERENCE
% ---------------------------------------------------------------

[A0,T] = readgeoraster(Filenames_PET(1));

[nRows,nCols] = size(A0);

fprintf('Original raster dimensions: %d x %d\n', ...
    nRows,nCols);

clear A0;


%% ---------------------------------------------------------------
% 8. CHECK REGION MASK
% ---------------------------------------------------------------

if ~isequal(size(RegionMask_PET),[nRows,nCols])

    error(['RegionMask_PET dimensions do not match ', ...
           'the MOD16 PET raster.']);

end

% Convert mask to logical
regionMask = RegionMask_PET ~= 0;


%% ---------------------------------------------------------------
% 9. PREALLOCATE QC INFORMATION
% ---------------------------------------------------------------

QC_Date = monthDates;

QC_Mean   = nan(nMonths,1);
QC_Median = nan(nMonths,1);
QC_Min    = nan(nMonths,1);
QC_Max    = nan(nMonths,1);
QC_P01    = nan(nMonths,1);
QC_P99    = nan(nMonths,1);

QC_ValidPixels = zeros(nMonths,1);
QC_UsedComposites = zeros(nMonths,1);


%% ---------------------------------------------------------------
% 10. PROCESS ONE MONTH AT A TIME
% ---------------------------------------------------------------

for im = 1:nMonths

    monthStart = monthDates(im);
    monthEnd = dateshift(monthStart,'end','month');

    fprintf('\n----------------------------------------\n');
    fprintf('Processing %s\n',datestr(monthStart,'yyyy-mm'));


    %% Identify all 8-day composites overlapping this month

    overlap = ...
        startDate <= monthEnd & ...
        endDate >= monthStart;

    idx = find(overlap);

    QC_UsedComposites(im) = numel(idx);

    fprintf('Overlapping composites: %d\n',numel(idx));


    %% Monthly accumulator

    MonthlyPET = zeros(nRows,nCols,'single');

    % Number of valid temporal contributions at each pixel
    ValidCount = zeros(nRows,nCols,'uint8');


    %% -----------------------------------------------------------
    % Process each overlapping composite
    % -----------------------------------------------------------

    for j = 1:numel(idx)

        k = idx(j);

        [RawPET,Rcurrent] = readgeoraster(Filenames_PET(k));

        RawPET = double(RawPET);


        %% -------------------------------------------------------
        % A. REMOVE MOD16 SPECIAL/FILL VALUES BEFORE SCALING
        % -------------------------------------------------------

        % MOD16 special values:
        % 32761 = unclassified
        % 32762 = urban/built-up
        % 32763 = wetland
        % 32764 = snow/ice
        % 32765 = barren/sparse vegetation
        % 32766 = water
        % 32767 = fill

        invalid = ...
            RawPET >= 32761 & RawPET <= 32767;

        RawPET(invalid) = NaN;


        %% Remove negative values for PET

        RawPET(RawPET < 0) = NaN;


        %% -------------------------------------------------------
        % B. APPLY MOD16 SCALE FACTOR
        % -------------------------------------------------------

        PET8 = RawPET * 0.1;

        % PET8 is now equivalent to mm over the composite period.


        %% -------------------------------------------------------
        % C. MASK OUTSIDE WEST AFRICA
        % -------------------------------------------------------

        PET8(~regionMask) = NaN;


        %% -------------------------------------------------------
        % D. NUMBER OF COMPOSITE DAYS FALLING IN CURRENT MONTH
        % -------------------------------------------------------

        overlapStart = max(startDate(k),monthStart);
        overlapEnd   = min(endDate(k),monthEnd);

        nOverlapDays = days(overlapEnd-overlapStart) + 1;

        totalCompositeDays = nDaysComposite(k);


        %% -------------------------------------------------------
        % E. TEMPORAL ALLOCATION
        %
        % PET8 is a total over the complete composite period.
        %
        % Fraction allocated to this calendar month:
        %
        % nOverlapDays / totalCompositeDays
        % -------------------------------------------------------

        fraction = nOverlapDays / totalCompositeDays;

        PETContribution = PET8 * fraction;


        %% -------------------------------------------------------
        % F. ADD ONLY VALID VALUES
        % -------------------------------------------------------

        valid = isfinite(PETContribution);

        MonthlyPET(valid) = ...
            MonthlyPET(valid) + single(PETContribution(valid));

        ValidCount(valid) = ValidCount(valid) + 1;


        clear RawPET PET8 PETContribution invalid valid

    end


    %% -----------------------------------------------------------
    % 11. REMOVE PIXELS WITHOUT VALID DATA
    % -----------------------------------------------------------

    MonthlyPET(ValidCount == 0) = NaN;

    MonthlyPET(~regionMask) = NaN;


    %% -----------------------------------------------------------
    % 12. QUALITY-CONTROL STATISTICS BEFORE RESAMPLING
    % -----------------------------------------------------------

    values = double(MonthlyPET(regionMask));
    values = values(isfinite(values));

    if isempty(values)

        warning('No valid PET pixels for %s.', ...
            datestr(monthStart,'yyyy-mm'));

    else

        QC_Mean(im)   = mean(values);
        QC_Median(im) = median(values);

        QC_Min(im) = min(values);
        QC_Max(im) = max(values);

        QC_P01(im) = prctile(values,1);
        QC_P99(im) = prctile(values,99);

        QC_ValidPixels(im) = numel(values);

        fprintf(['Mean = %.2f | Median = %.2f | ', ...
                 'P01 = %.2f | P99 = %.2f | ', ...
                 'Min = %.2f | Max = %.2f\n'], ...
                 QC_Mean(im), ...
                 QC_Median(im), ...
                 QC_P01(im), ...
                 QC_P99(im), ...
                 QC_Min(im), ...
                 QC_Max(im));

    end


    %% -----------------------------------------------------------
    % 13. SAVE MONTHLY ORIGINAL-GRID PET
    % -----------------------------------------------------------

    dateString = string(datestr(monthStart,'yyyy.mm'));

    Outname = "PET.Monthly.Corrected." + dateString + ".tif";

    outputFile = fullfile(TemporaryOutDir,Outname);

    geotiffwrite(outputFile,MonthlyPET,T);


    %% -----------------------------------------------------------
    % 14. RESAMPLE TO EXISTING TARGET GRID
    %
    % This follows your original GMT grdtrack workflow.
    %
    % IMPORTANT:
    % Resampling is done AFTER fill-value masking and monthly
    % temporal aggregation.
    % -----------------------------------------------------------
    [m, n] = size(NDVI_LonGrid);
    G = gmt(char("grdtrack -N -G" + outputFile),[NDVI_LonGrid(:),NDVI_LatGrid(:)]);

    ResampledPET = G.data(:,3);

    ResampledPET = reshape(ResampledPET,m,n);


    %% Remove non-finite output

    ResampledPET(~isfinite(ResampledPET)) = NaN;


    %% Save resampled monthly PET

    Outname2 =  "PET.Monthly.Resampled.Corrected." + dateString + ".tif";

    outputFile2 = fullfile(OutDir,Outname2);

    geotiffwrite(outputFile2,flipud(ResampledPET),R);


    %% Free memory
    

    clear MonthlyPET ValidCount values ...
          G ResampledPET

end


%% ---------------------------------------------------------------
% 15. EXPORT PROCESSING QC TO EXCEL
% ---------------------------------------------------------------

QC_Table = table( ...
    QC_Date, ...
    QC_UsedComposites, ...
    QC_ValidPixels, ...
    QC_Mean, ...
    QC_Median, ...
    QC_Min, ...
    QC_Max, ...
    QC_P01, ...
    QC_P99, ...
    'VariableNames', ...
    {'Date','N_Composites','ValidPixels', ...
     'MeanPET_mm_month','MedianPET_mm_month', ...
     'MinPET_mm_month','MaxPET_mm_month', ...
     'P01_mm_month','P99_mm_month'});

QCfile = fullfile(OutDir,'PET_Monthly_Processing_QC.xlsx');

writetable(QC_Table,QCfile);

fprintf('\n============================================\n');
fprintf('PET PROCESSING COMPLETED\n');
fprintf('============================================\n');
fprintf('QC file:\n%s\n',QCfile);
%%
%% ================================================================
% 16. NORMALIZE RESAMPLED MONTHLY PET
% ================================================================
%
% Pixel-wise min-max normalization over the complete study period:
%
%              PET(t) - PETmin
%   NPET(t) = -----------------
%              PETmax - PETmin
%
% PETmin and PETmax are calculated independently for each pixel
% using all resampled monthly PET images from 2001-2022.
%
% Output range:
%   NPET = 0  -> temporal minimum PET at that pixel
%   NPET = 1  -> temporal maximum PET at that pixel
%
% ================================================================


%% ---------------------------------------------------------------
% 16.1. DEFINE OUTPUT FOLDER FOR NORMALIZED PET
% ---------------------------------------------------------------

NormalizedOutDir = ...
"G:\PHD_Folders\PhD_Data\West_Africa_data\All_data_for_For_Thesis\Richard_Results\PET.Monthly.Normalized.Corrected";

if ~exist(NormalizedOutDir,'dir')
    mkdir(NormalizedOutDir);
end

fprintf('\n============================================\n');
fprintf('STARTING PET NORMALIZATION\n');
fprintf('============================================\n');


%% ---------------------------------------------------------------
% 16.2. LIST ALL RESAMPLED MONTHLY PET FILES
% ---------------------------------------------------------------

DDD_Resampled = dir(fullfile(OutDir, ...
    'PET.Monthly.Resampled.Corrected.*.tif'));

if isempty(DDD_Resampled)
    error('No resampled PET GeoTIFF files were found.');
end

Names_Resampled = {DDD_Resampled.name}';
Names_Resampled = string(Names_Resampled);

Folder_Resampled = {DDD_Resampled.folder}';
Folder_Resampled = string(Folder_Resampled);

Files_Resampled = fullfile( ...
    Folder_Resampled, ...
    Names_Resampled);

fprintf('Found %d resampled monthly PET files.\n', ...
    numel(Files_Resampled));


%% ---------------------------------------------------------------
% 16.3. CHECK EXPECTED NUMBER OF MONTHS
% ---------------------------------------------------------------

if numel(Files_Resampled) ~= 264

    warning(['Expected 264 monthly PET files for 2001-2022, ', ...
             'but found %d files.'], ...
             numel(Files_Resampled));

end


%% ---------------------------------------------------------------
% 16.4. READ FIRST RESAMPLED PET IMAGE
% ---------------------------------------------------------------

[FirstPET,R_PET] = readgeoraster(Files_Resampled(1));

FirstPET = double(FirstPET);

[nRowsPET,nColsPET] = size(FirstPET);

fprintf('Resampled PET raster dimensions: %d x %d\n', ...
    nRowsPET,nColsPET);


%% ---------------------------------------------------------------
% 16.5. INITIALIZE PIXEL-WISE TEMPORAL MINIMUM AND MAXIMUM
% ---------------------------------------------------------------

% Initialize with NaN.
PET_Min = nan(nRowsPET,nColsPET);
PET_Max = nan(nRowsPET,nColsPET);


%% ---------------------------------------------------------------
% 16.6. FIRST PASS:
% CALCULATE TEMPORAL MINIMUM AND MAXIMUM FOR EACH PIXEL
% ---------------------------------------------------------------

fprintf('\nCalculating pixel-wise temporal PET minimum and maximum...\n');

for i = 1:numel(Files_Resampled)

    PET = readgeoraster(Files_Resampled(i));

    PET = double(PET);

    % Make sure non-finite values remain NaN
    PET(~isfinite(PET)) = NaN;


    %% Identify valid PET pixels

    valid = isfinite(PET);


    %% Initialize PET_Min and PET_Max where necessary

    initializeMin = valid & isnan(PET_Min);
    initializeMax = valid & isnan(PET_Max);

    PET_Min(initializeMin) = PET(initializeMin);
    PET_Max(initializeMax) = PET(initializeMax);


    %% Update temporal minimum

    updateMin = valid & ...
                (isnan(PET_Min) | PET < PET_Min);

    PET_Min(updateMin) = PET(updateMin);


    %% Update temporal maximum

    updateMax = valid & ...
                (isnan(PET_Max) | PET > PET_Max);

    PET_Max(updateMax) = PET(updateMax);


    fprintf('Min/Max pass: %d of %d\n', ...
        i,numel(Files_Resampled));

end


%% ---------------------------------------------------------------
% 16.7. CALCULATE TEMPORAL PET RANGE
% ---------------------------------------------------------------

PET_Range = PET_Max - PET_Min;

% Pixels for which PETmax = PETmin cannot be normalized
ZeroRange = PET_Range == 0;

PET_Range(ZeroRange) = NaN;

fprintf('\nPixels with PETmax = PETmin: %d\n', ...
    sum(ZeroRange(:)));


%% ---------------------------------------------------------------
% 16.8. SECOND PASS:
% NORMALIZE EACH MONTH AND SAVE AS GEOTIFF
% ---------------------------------------------------------------

fprintf('\nNormalizing monthly PET images...\n');

for i = 1:numel(Files_Resampled)

    PET = readgeoraster(Files_Resampled(i));

    PET = double(PET);

    PET(~isfinite(PET)) = NaN;


    %% -----------------------------------------------------------
    % PIXEL-WISE MIN-MAX NORMALIZATION
    % -----------------------------------------------------------

    NPET = ...
        (PET - PET_Min) ./ PET_Range;


    %% -----------------------------------------------------------
    % Remove non-finite values
    % -----------------------------------------------------------

    NPET(~isfinite(NPET)) = NaN;


    %% -----------------------------------------------------------
    % Numerical protection
    %
    % The theoretical range is [0,1].
    % These two lines only protect against very small floating-point
    % rounding errors.
    % -----------------------------------------------------------

    NPET(NPET < 0) = 0;
    NPET(NPET > 1) = 1;


    %% -----------------------------------------------------------
    % EXTRACT DATE FROM ORIGINAL RESAMPLED FILENAME
    %
    % Input example:
    % PET.Monthly.Resampled.Corrected.2001.01.tif
    %
    % Output:
    % NPET.Monthly.2001.01.tif
    % -----------------------------------------------------------

    CurrentName = Names_Resampled(i);

    DatePart = erase( ...
        CurrentName, ...
        "PET.Monthly.Resampled.Corrected.");

    DatePart = erase(DatePart,'.tif');


    %% -----------------------------------------------------------
    % OUTPUT FILENAME
    % -----------------------------------------------------------

    NormalizedName = ...
        "NPET.Monthly." + DatePart + ".tif";

    NormalizedFile = ...
        fullfile(NormalizedOutDir,NormalizedName);


    %% -----------------------------------------------------------
    % SAVE NORMALIZED PET
    % -----------------------------------------------------------

    geotiffwrite( ...
        NormalizedFile, ...
        NPET, ...
        R_PET);


    fprintf('Normalized PET: %d of %d -> %s\n', ...
        i,numel(Files_Resampled),NormalizedName);

end


%% ---------------------------------------------------------------
% 16.9. SAVE TEMPORAL MINIMUM AND MAXIMUM PET RASTERS
% ---------------------------------------------------------------
%
% These files are useful for reproducibility because they document
% the pixel-specific bounds used in the normalization.
% ---------------------------------------------------------------

MinFile = fullfile( ...
    NormalizedOutDir, ...
    'PET_Temporal_Min_2001_2022.tif');

MaxFile = fullfile( ...
    NormalizedOutDir, ...
    'PET_Temporal_Max_2001_2022.tif');


geotiffwrite( ...
    MinFile, ...
    PET_Min, ...
    R_PET);

geotiffwrite( ...
    MaxFile, ...
    PET_Max, ...
    R_PET);


%% ---------------------------------------------------------------
% 16.10. FINAL CHECK
% ---------------------------------------------------------------

DDD_NPET = dir(fullfile( ...
    NormalizedOutDir, ...
    'NPET.Monthly.*.tif'));

fprintf('\n============================================\n');
fprintf('PET NORMALIZATION COMPLETED\n');
fprintf('============================================\n');

fprintf('Normalized PET files created: %d\n', ...
    numel(DDD_NPET));

fprintf('Output folder:\n%s\n', ...
    NormalizedOutDir);

fprintf('\nTemporal minimum PET:\n%s\n', ...
    MinFile);

fprintf('\nTemporal maximum PET:\n%s\n', ...
    MaxFile);

fprintf('============================================\n');