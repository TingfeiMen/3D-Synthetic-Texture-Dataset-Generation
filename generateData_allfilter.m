%GENERATEDATA_ALLFILTER Generate a 3-D synthetic texture dataset and its Haralick features.
%
%   This script builds a labelled dataset of synthetic 3-D textures. For each
%   wavelet filter in WAVELETLIST it draws NUMSAMPLES white Gaussian noise
%   volumes, decomposes each volume with WAVEDEC3, and then enumerates every
%   way of zeroing out N of the resulting subbands. Each such combination is
%   reconstructed with WAVEREC3 and treated as one texture class: the
%   surviving subbands determine the spectral signature of the volume.
%
%   Every reconstructed volume is quantized to 8 bits, saved as a .mat file,
%   and summarized by 13 Haralick features computed from its grey-level
%   co-occurrence matrix. Features are appended to a per-sample CSV file.
%
%   Layout produced under OUTPUTROOT
%     TestData_AllFilters_Var<variance>_<numSamples>/
%       Wname_<wavelet>/
%         FeatureData/
%           ReconstructedData_<sampleIdx>.csv     features, one row per class
%         ReconstructedData_<sampleIdx>/
%           Sample_<wavelet>_Var<variance>_<zeroedSubbands>.mat
%
%   The script has two sections. Section 1 generates and writes the dataset.
%   Section 2 is a standalone benchmark that repeats the same pipeline with a
%   smaller sample count and writes nothing to disk; it reports how the
%   Haralick cost splits between GLCM construction and feature extraction.
%   Run the sections independently with Ctrl+Enter.
%
%   Requirements
%     MATLAB with Wavelet Toolbox (WAVEDEC3, WAVEREC3).
%     Local functions of this project: NORMALIZEANDSCALE, CALCULATE3DGLCM,
%     COMPUTEHARALICK, all on the MATLAB path.
%
%   Configuration
%     Edit the parameter block below. OUTPUTROOT is the only path that needs
%     to change when moving between machines.
%
%   WARNING - disk usage
%     The defaults are large. With numSamples = 70, six wavelets and
%     n = 6 zeroed subbands, NCHOOSEK(numCoeffs, 6) runs to tens of thousands
%     of classes per sample, and each 64x64x64 uint8 volume is 256 KB. A full
%     run at these settings produces on the order of hundreds of gigabytes.
%     Reduce NUMSAMPLES or use the Section 2 benchmark before committing to a
%     full run, and point OUTPUTROOT at a drive with room to spare rather
%     than at a cloud-synced folder.
%
%   Reproducibility
%     RNG("default") is called once at the start of each section, so a given
%     parameter set reproduces the same volumes. Changing NUMSAMPLES,
%     DIMENSIONS or the order of WAVELETLIST changes the draw sequence and
%     therefore the data.
%
%   Note on GLCM normalization
%     CALCULATE3DGLCM contains a known normalization defect that is retained
%     deliberately for reproducibility. See the WARNING in that function
%     before reusing the feature values in new work.
%
%   See also NORMALIZEANDSCALE, CALCULATE3DGLCM, COMPUTEHARALICK,
%            WAVEDEC3, WAVEREC3, NCHOOSEK.

%   Authors:     Tingfei Men, Ahmed Ashraf, Sherif S. Sherif
%   Affiliation: Department of Electrical and Computer Engineering,
%                University of Manitoba, Winnipeg, Manitoba, Canada
%   Contact:     Tingfei Men <tingfeimen@gmail.com>
%
%   Citation:    If you use this code in a scientific publication, please cite
%                T. Men, A. Ashraf, and S. S. Sherif, "Generation of a large
%                three-dimensional texture image dataset for machine and deep
%                learning applications," accepted to the International Workshop
%                on New Approaches to Multivariate Statistical Analysis
%                (NAMSP 2026).
%   Project:     3D Synthetic Texture Dataset Generation
%   Created:     2025-10-08
%   Last update: 2026-08-17
%   License:     MIT. See the LICENSE file in the repository root.

%% Section 1 - Generate the dataset and extract Haralick features

clear;
clc;

% ----- Output location -------------------------------------------------
% Destination for the generated dataset. Edit this for your machine; it is
% created if it does not exist. Prefer a local drive with plenty of free
% space over a cloud-synced folder.
outputRoot = fullfile(pwd, 'data');

% ----- Sampling parameters ---------------------------------------------
numSamples = 70;            % Noise volumes drawn per wavelet filter
n = 6;                      % Number of wavelet subbands zeroed per class
mean_value = 0;             % Mean of the white Gaussian noise
variance = 10;              % Variance of the white Gaussian noise
stddev = sqrt(variance);
dimensions = [64, 64, 64];  % Size of each generated volume
rng("default");             % Fixed seed for reproducibility

% ----- Wavelet decomposition parameters --------------------------------
waveletList = {'db1', 'db2', 'db3', 'coif1', 'sym4', 'fk4'};
levels = 2;                 % Decomposition depth passed to WAVEDEC3

% ----- GLCM settings ---------------------------------------------------
% The 13 unique directions of a 3-D neighbourhood at unit distance.
directions = [
    1, 0, 0;  0, 1, 0;  0, 0, 1;
    1, 1, 0;  1, -1, 0;
    1, 0, 1;  1, 0, -1;
    0, 1, 1;  0, 1, -1;
    1, 1, 1;  1, -1, 1;
    1, 1, -1; 1, -1, -1
];
grayLevels = 256;           % Matches the uint8 output of normalizeAndScale
d = 1;                      % Index of the direction whose GLCM is featurized

% Root folder for all filters
mainRootDir = fullfile(outputRoot, ...
    sprintf('TestData_AllFilters_Var%d_%d', variance, numSamples));
if ~exist(mainRootDir, 'dir')
    mkdir(mainRootDir);
end

totalHaralickTime = 0;
tic

for w = 1:length(waveletList)
    wname = waveletList{w};

    % Folder setup for each wavelet filter
    mainDataDir = fullfile(mainRootDir, sprintf('Wname_%s', wname));
    featureDataDir = fullfile(mainDataDir, 'FeatureData');
    if ~exist(mainDataDir, 'dir')
        mkdir(mainDataDir);
    end
    if ~exist(featureDataDir, 'dir')
        mkdir(featureDataDir);
    end

    for sampleIdx = 1:numSamples
        outputDir = fullfile(mainDataDir, sprintf('ReconstructedData_%d', sampleIdx));
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        % One CSV per sample. The column order below must match the output
        % vector of computeHaralick; see that function's help text.
        sampleCsvFilename = fullfile(featureDataDir, sprintf('ReconstructedData_%d.csv', sampleIdx));
        header = ["ReconstructedCube", "Wavelet", "Energy", "Correlation", "Inertia", "Entropy", ...
                  "Inverse Diff", "Sum Avg", "Sum Var", "Sum Entropy", ...
                  "Diff Avg", "Diff Var", "Diff Entropy", "Info Corr 1", "Info Corr 2"];
        writecell(cellstr(header), sampleCsvFilename, 'WriteMode', 'overwrite');

        % Generate random sample
        currentSample = mean_value + stddev * randn(dimensions);

        % Wavelet decomposition
        WD = wavedec3(currentSample, levels, wname);
        numCoeffs = length(WD.dec);

        % Each row lists the n subbands zeroed out for one texture class
        zeroOutCombinations = nchoosek(1:numCoeffs, n);
        numClasses = size(zeroOutCombinations, 1);

        coeffMatrix = cell(numCoeffs, 1);
        for i = 1:numCoeffs
            coeffMatrix{i} = WD.dec{i};
        end

        % Loop over zeroing combinations
        for classIdx = 1:numClasses
            tempCoeffMatrix = coeffMatrix;
            zeroedOut = zeroOutCombinations(classIdx, :);

            for j = 1:n
                tempCoeffMatrix{zeroedOut(j)} = zeros(size(tempCoeffMatrix{zeroedOut(j)}));
            end

            % Reconstruct
            WD_mod = WD;
            WD_mod.dec = tempCoeffMatrix;
            reconstructedSample = waverec3(WD_mod);
            sampleData = normalizeAndScale(reconstructedSample);

            % Save .mat file. The zeroed subband indices are encoded in the
            % filename and serve as the class label.
            zeroOutString = sprintf('_%d', zeroedOut);
            matFilename = fullfile(outputDir, sprintf('Sample_%s_Var%d%s.mat', wname, variance, zeroOutString));
            save(matFilename, 'sampleData');

            % Haralick feature extraction
            haralickStart = tic;
            glcm = calculate3DGLCM(sampleData, directions, grayLevels);
            F = computeHaralick(glcm(:, :, d));
            haralickElapsed = toc(haralickStart);
            totalHaralickTime = totalHaralickTime + haralickElapsed;

            % Write to CSV
            reconstructedFilename = sprintf('Sample_%s_Var%d%s.mat', wname, variance, zeroOutString);
            row = [{reconstructedFilename}, {wname}, num2cell(F)];
            writecell(row, sampleCsvFilename, 'WriteMode', 'append');
        end
    end
end

toc
fprintf('Total Haralick computation time: %.2f seconds\n', totalHaralickTime);

%% Section 2 - Benchmark only, writes nothing to disk
%  Repeats the pipeline with a reduced sample count and separates the cost of
%  building the GLCM from the cost of extracting features. Use this to size a
%  full run before starting one.

clear;
clc;

% ----- Sampling parameters (reduced for benchmarking) ------------------
numSamples = 7;
n = 6;
mean_value = 0;
variance = 10;
stddev = sqrt(variance);
dimensions = [64, 64, 64];
rng("default");

% ----- Wavelet decomposition parameters --------------------------------
waveletList = {'db1', 'db2', 'db3', 'coif1', 'sym4', 'fk4'};
levels = 2;

% ----- GLCM settings ---------------------------------------------------
directions = [
    1, 0, 0;  0, 1, 0;  0, 0, 1;
    1, 1, 0;  1, -1, 0;
    1, 0, 1;  1, 0, -1;
    0, 1, 1;  0, 1, -1;
    1, 1, 1;  1, -1, 1;
    1, 1, -1; 1, -1, -1
];
grayLevels = 256;
d = 1;

totalHaralickTime = 0;
totalGLCMTime     = 0;
totalFeatureTime  = 0;
totalSamples      = 0;

disp('Starting Haralick computation...');

for w = 1:length(waveletList)
    wname = waveletList{w};

    for sampleIdx = 1:numSamples
        currentSample = mean_value + stddev * randn(dimensions);

        WD = wavedec3(currentSample, levels, wname);
        numCoeffs = length(WD.dec);
        zeroOutCombinations = nchoosek(1:numCoeffs, n);
        numClasses = size(zeroOutCombinations, 1);

        coeffMatrix = cell(numCoeffs, 1);
        for i = 1:numCoeffs
            coeffMatrix{i} = WD.dec{i};
        end

        for classIdx = 1:numClasses
            tempCoeffMatrix = coeffMatrix;
            zeroedOut = zeroOutCombinations(classIdx, :);
            for j = 1:n
                tempCoeffMatrix{zeroedOut(j)} = zeros(size(tempCoeffMatrix{zeroedOut(j)}));
            end

            WD_mod = WD;
            WD_mod.dec = tempCoeffMatrix;
            reconstructedSample = waverec3(WD_mod);
            sampleData = normalizeAndScale(reconstructedSample);

            % --- Pure GLCM timing ---
            glcmStart = tic;
            glcm = calculate3DGLCM(sampleData, directions, grayLevels);
            glcmElapsed = toc(glcmStart);
            totalGLCMTime = totalGLCMTime + glcmElapsed;

            % --- Pure feature extraction timing ---
            featureStart = tic;
            F = computeHaralick(glcm(:, :, d));
            featureElapsed = toc(featureStart);
            totalFeatureTime = totalFeatureTime + featureElapsed;

            totalHaralickTime = totalHaralickTime + glcmElapsed + featureElapsed;
            totalSamples = totalSamples + 1;
        end
    end
end

fprintf('\n========== Haralick Timing Summary ==========\n');
fprintf('Total Pure Haralick Time:         %.4f seconds\n', totalHaralickTime);
fprintf('  |- GLCM Construction:           %.4f seconds\n', totalGLCMTime);
fprintf('  |- Feature Extraction:          %.4f seconds\n', totalFeatureTime);
fprintf('Wavelets Processed:               %d\n',           length(waveletList));
fprintf('Total Cubes Processed:            %d\n',           totalSamples);
fprintf('Avg Haralick Time per Cube:       %.4f ms\n',      totalHaralickTime / totalSamples * 1000);
fprintf('==============================================\n');
