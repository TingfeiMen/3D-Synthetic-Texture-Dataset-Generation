function GLCMs = calculate3DGLCM(matrix, directions, grayLevels)
%CALCULATE3DGLCM Grey-level co-occurrence matrices of a 3-D volume.
%
%   GLCMS = CALCULATE3DGLCM(MATRIX, DIRECTIONS, GRAYLEVELS) computes one
%   grey-level co-occurrence matrix (GLCM) per direction for the 3-D
%   integer-valued volume MATRIX. For each direction, every voxel is paired
%   with the neighbour reached by that offset, and the corresponding entry
%   of the co-occurrence matrix is incremented. Voxel pairs whose neighbour
%   falls outside the volume are discarded.
%
%   Inputs
%     matrix      - 3-D array of integer grey levels in [0, GRAYLEVELS-1],
%                   e.g. the uint8 output of NORMALIZEANDSCALE.
%     directions  - N-by-3 array. Each row is a voxel offset [dx dy dz]
%                   applied along the first, second and third dimension
%                   respectively, e.g. [1 0 0] or [1 -1 1]. The 13 unique
%                   directions of a 3-D neighbourhood at distance 1 are the
%                   usual choice.
%     grayLevels  - Number of grey levels, e.g. 256 for uint8 data. Sets the
%                   size of each co-occurrence matrix.
%
%   Outputs
%     GLCMs       - GRAYLEVELS-by-GRAYLEVELS-by-N array. GLCMS(:,:,K) holds
%                   the co-occurrence matrix for DIRECTIONS(K,:).
%
%   Example
%     directions = [1 0 0; 0 1 0; 0 0 1];
%     vol        = normalizeAndScale(randn(64, 64, 64));
%     glcm       = calculate3DGLCM(vol, directions, 256);
%     F          = computeHaralick(glcm(:, :, 1));
%
%   Notes
%     - The co-occurrence counts are directed, not symmetric: the pair
%       (a, b) along [1 0 0] is counted, but the reverse pair (b, a) is not
%       added. Symmetric GLCMs would require adding the transpose, or
%       including both a direction and its negation in DIRECTIONS.
%     - Grey level g is stored at index g+1, since MATLAB indexes from 1.
%     - The implementation uses explicit nested loops for clarity. Cost grows
%       as NUMEL(MATRIX) * N, so large volumes or many directions are slow.
%
%   WARNING - known normalization defect, retained deliberately
%     The normalization statement at the end of the direction loop divides
%     the ENTIRE GLCMS array, and it executes once per direction. As a
%     result GLCMS(:,:,1) is divided N times, GLCMS(:,:,2) is divided N-1
%     times, and so on, so the returned matrices are NOT each normalized to
%     unit sum and are not comparable across directions. Only the last
%     direction is scaled exactly once.
%
%     This behaviour is preserved intentionally so that results remain
%     reproducible against datasets already generated with this code. A
%     corrected version would move the normalization out of the loop and
%     normalize each direction separately, for example:
%
%         for k = 1:numDirections
%             GLCMs(:,:,k) = GLCMs(:,:,k) / sum(sum(GLCMs(:,:,k)));
%         end
%
%     Downstream Haralick features assume a GLCM that sums to 1, so anyone
%     reusing this function for new work should apply the correction above
%     and regenerate their features.
%
%   See also NORMALIZEANDSCALE, COMPUTEHARALICK, GRAYCOMATRIX.

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
%   Created:     2025-01-13
%   Last update: 2026-08-17
%   License:     MIT. See the LICENSE file in the repository root.

    % Get the size of the matrix
    [rows, cols, depths] = size(matrix);

    % Number of directions
    numDirections = size(directions, 1);

    % Initialize GLCMs
    GLCMs = zeros(grayLevels, grayLevels, numDirections);

    % Loop through each direction
    for dirIdx = 1:numDirections
        % Get the current direction (offset)
        offset = directions(dirIdx, :);
        [dx, dy, dz] = deal(offset(1), offset(2), offset(3));

        % Initialize the GLCM for this direction
        glcm = zeros(grayLevels, grayLevels);

        % Iterate over the matrix
        for z = 1:depths
            for y = 1:cols
                for x = 1:rows
                    % Calculate neighbor indices
                    nbrX = x + dx;
                    nbrY = y + dy;
                    nbrZ = z + dz;

                    % Ensure neighbor indices are within bounds
                    if nbrX > 0 && nbrX <= rows && nbrY > 0 && nbrY <= cols && nbrZ > 0 && nbrZ <= depths
                        % Get the reference pixel and the neighbor pixel
                        refValue = matrix(x, y, z);
                        nbrValue = matrix(nbrX, nbrY, nbrZ);

                        % Update the GLCM
                        glcm(refValue + 1, nbrValue + 1) = glcm(refValue + 1, nbrValue + 1) + 1;
                    end
                end
            end
        end

        % Store the GLCM for this direction
        GLCMs(:, :, dirIdx) = glcm;
        % summed up to get direction-invariant GLCM
        % GLCMs = sum(GLCMs,3);

        % See the WARNING in the header: this normalizes the whole array and
        % runs once per direction. Kept unchanged for reproducibility.
        GLCMs = GLCMs./sum(sum(GLCMs));
    end
end
