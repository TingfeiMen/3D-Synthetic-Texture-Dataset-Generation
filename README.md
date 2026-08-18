# 3D Synthetic Texture Dataset Generation

MATLAB code for generating a large dataset of labelled three-dimensional
synthetic texture volumes, and for extracting Haralick texture features from
them.

## Abstract

Machine and deep learning models for volumetric texture analysis are limited by
the scarcity of large, labelled 3-D texture datasets. This code generates such
a dataset synthetically and at scale.

The generator draws white Gaussian noise volumes and decomposes each one with a
three-dimensional discrete wavelet transform. It then enumerates every way of
zeroing out a fixed number of the resulting subbands and reconstructs a volume
from each surviving subset. Because the retained subbands determine the
spectral content of the reconstruction, each combination yields a distinct and
exactly known texture class, and the combination itself serves as the class
label. Repeating this over several wavelet families and many noise
realizations produces a large dataset with controlled, reproducible ground
truth.

Each reconstructed volume is quantized to 8 bits and summarized by 13 Haralick
texture features computed from its three-dimensional grey-level co-occurrence
matrix (GLCM).

## Requirements

- MATLAB (developed and tested on R2024a and later)
- Wavelet Toolbox, for `wavedec3` and `waverec3`

## Repository contents

| File | Description |
| --- | --- |
| `generateData_allfilter.m` | Main script. Generates the dataset and writes the feature CSV files. |
| `calculate3DGLCM.m` | Grey-level co-occurrence matrices of a 3-D volume, one per direction. |
| `computeHaralick.m` | 13 Haralick texture features from a 2-D co-occurrence matrix. |
| `normalizeAndScale.m` | Min-max normalization of an array to the 8-bit range `[0, 255]`. |

## Example

Clone the repository, add it to the MATLAB path, and open the main script:

```matlab
addpath(genpath('3D-Synthetic-Texture-Dataset-Generation'));
edit generateData_allfilter
```

Set the output location near the top of Section 1, then run that section:

```matlab
outputRoot = fullfile(pwd, 'data');   % edit for your machine
```

The individual functions can also be used on their own:

```matlab
% One 64x64x64 volume, quantized and featurized along the [1 0 0] direction
directions = [1 0 0; 0 1 0; 0 0 1];
vol        = normalizeAndScale(randn(64, 64, 64));
glcm       = calculate3DGLCM(vol, directions, 256);
F          = computeHaralick(glcm(:, :, 1));
```

Section 2 of `generateData_allfilter.m` is a self-contained benchmark. It runs
the same pipeline with a reduced sample count, writes nothing to disk, and
reports how the cost splits between building the GLCM and extracting features.
Run it first to estimate the cost of a full run.

## Function calls

```matlab
matrixScaled = normalizeAndScale(matrix)
GLCMs        = calculate3DGLCM(matrix, directions, grayLevels)
F            = computeHaralick(GLCM)
```

### `normalizeAndScale`

| Argument | Description |
| --- | --- |
| `matrix` | Real-valued 2-D or 3-D numeric array. |
| **Returns** `matrixScaled` | `uint8` array of the same size, values in `[0, 255]`. |

### `calculate3DGLCM`

| Argument | Description |
| --- | --- |
| `matrix` | 3-D array of integer grey levels in `[0, grayLevels-1]`. |
| `directions` | `N`-by-3 array of voxel offsets `[dx dy dz]`. The 13 unique directions of a 3-D unit neighbourhood are the usual choice. |
| `grayLevels` | Number of grey levels, e.g. `256` for `uint8` data. |
| **Returns** `GLCMs` | `grayLevels`-by-`grayLevels`-by-`N` array; slice `k` corresponds to `directions(k,:)`. |

Co-occurrence counts are directed, not symmetric: the pair `(a, b)` along
`[1 0 0]` is counted, but the reverse pair is not added. Voxel pairs whose
neighbour falls outside the volume are discarded.

### `computeHaralick`

| Argument | Description |
| --- | --- |
| `GLCM` | `G`-by-`G` normalized co-occurrence matrix, e.g. one slice of `calculate3DGLCM`. |
| **Returns** `F` | 1-by-13 row vector of features. |

Features are returned in this order, which matches the CSV header written by
`generateData_allfilter.m`:

| # | Feature | # | Feature |
| --- | --- | --- | --- |
| 1 | Energy (angular second moment) | 8 | Sum entropy |
| 2 | Correlation | 9 | Difference average |
| 3 | Inertia (contrast) | 10 | Difference variance |
| 4 | Entropy | 11 | Difference entropy |
| 5 | Inverse difference moment | 12 | Information measure of correlation 1 |
| 6 | Sum average | 13 | Information measure of correlation 2 |
| 7 | Sum variance | | |

Only 13 features are computed. Haralick's original set has 14; the maximal
correlation coefficient, which requires an eigenvalue decomposition, is
omitted.

## Output layout

```
<outputRoot>/
  TestData_AllFilters_Var<variance>_<numSamples>/
    Wname_<wavelet>/
      FeatureData/
        ReconstructedData_<sampleIdx>.csv      one row per texture class
      ReconstructedData_<sampleIdx>/
        Sample_<wavelet>_Var<variance>_<zeroedSubbands>.mat
```

The zeroed subband indices are encoded in each `.mat` filename and act as the
class label.

## Dataset examples

The `Dataset examples/` folder holds a small curated sample of the generated
dataset, so the data format can be inspected without regenerating anything.
It is a subset, not the full dataset.

```
Dataset examples/
  Wname_<wavelet>/                       db1, db2, db3, coif1, sym4, fk4
    FeatureData/
      ReconstructedData_1.csv            complete feature tables, all classes
      ...
      ReconstructedData_7.csv
    ReconstructedData/
      ReconstructedData_1/
        Sample_<wavelet>_Var10_<zeroedSubbands>.mat    10 sample volumes
```

For each of the six wavelet families:

- **All seven feature CSVs are included and complete.** Each file corresponds
  to one noise realization and holds one row per texture class, with the 13
  Haralick features in the order listed above. Nothing is truncated.
- **Ten sample volumes** are included from the first noise realization, taken
  at even intervals through the sorted class list so they span a range of
  subband combinations rather than clustering together.

Each `.mat` file stores a single variable `sampleData`, a 64x64x64 `uint8`
array. The zeroed subband indices are encoded in the filename and act as the
class label. To load one:

```matlab
S = load('Dataset examples/Wname_coif1/ReconstructedData/ReconstructedData_1/Sample_coif1_Var10_1_10_11_12_13_14.mat');
vol = S.sampleData;          % 64x64x64 uint8
volshow(vol);
```

## Reproducibility

`rng("default")` is called once at the start of each section, so a given
parameter set reproduces the same volumes. Changing `numSamples`, `dimensions`
or the order of `waveletList` changes the random draw sequence and therefore
the generated data.

## Citation

If you use this code in a scientific publication, please cite:

> T. Men, A. Ashraf, and S. S. Sherif, "Generation of a large three-dimensional
> texture image dataset for machine and deep learning applications," accepted to
> the International Workshop on New Approaches to Multivariate Statistical
> Analysis (NAMSP 2026).

## References

R. M. Haralick, K. Shanmugam and I. Dinstein, "Textural Features for Image
Classification," *IEEE Transactions on Systems, Man, and Cybernetics*,
vol. SMC-3, no. 6, pp. 610-621, 1973.

## Authors

Tingfei Men, Ahmed Ashraf, and Sherif S. Sherif  
Department of Electrical and Computer Engineering,  
University of Manitoba, Winnipeg, Manitoba, Canada

Corresponding author: Tingfei Men <tingfeimen@gmail.com>

## License

Released under the MIT License. See [LICENSE](LICENSE).
