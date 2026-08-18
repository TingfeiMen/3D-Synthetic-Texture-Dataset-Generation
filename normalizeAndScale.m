function matrixScaled = normalizeAndScale(matrix)
%NORMALIZEANDSCALE Linearly rescale an array to the 8-bit range [0, 255].
%
%   MATRIXSCALED = NORMALIZEANDSCALE(MATRIX) min-max normalizes MATRIX to
%   [0, 1] and then scales it to [0, 255], returning a uint8 array. This is
%   the quantization step applied to reconstructed wavelet volumes before
%   grey-level co-occurrence analysis, which requires integer grey levels.
%
%   Inputs
%     matrix        - Real-valued 2-D or 3-D numeric array. Typically a
%                     volume produced by RANDN or by WAVEREC3.
%
%   Outputs
%     matrixScaled  - uint8 array, same size as MATRIX, values in [0, 255].
%
%   Example
%     vol    = randn(64, 64, 64);
%     volU8  = normalizeAndScale(vol);
%     glcm   = calculate3DGLCM(volU8, [1 0 0], 256);
%
%   Notes
%     - The mapping is per-array, not per-slice: the global minimum and
%       maximum of MATRIX define the range, so the scaling is consistent
%       across the whole volume.
%     - Conversion to uint8 rounds to the nearest integer, so the operation
%       is lossy and not exactly invertible.
%     - If MATRIX is constant, maxValue equals minValue and the division
%       yields NaN; uint8(NaN) is 0, so a constant input maps to all zeros.
%
%   See also CALCULATE3DGLCM, COMPUTEHARALICK, UINT8, WAVEREC3.

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
%   Created:     2024-12-03
%   Last update: 2026-08-17
%   License:     MIT. See the LICENSE file in the repository root.

    % Step 1: Find the minimum and maximum values of the matrix
    minValue = min(matrix(:)); % Minimum value in the matrix
    maxValue = max(matrix(:)); % Maximum value in the matrix

    % Step 2: Normalize the matrix to the range [0, 1]
    matrixNormalized = (matrix - minValue) / (maxValue - minValue);

    % Step 3: Scale the normalized matrix to the range [0, 255] and convert to integer
    matrixScaled = uint8(matrixNormalized * 255);

end
