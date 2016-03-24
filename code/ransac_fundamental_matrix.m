% RANSAC Stencil Code
% CS 4495 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Find the best fundamental matrix using RANSAC on potentially matching
% points

% 'matches_a' and 'matches_b' are the Nx2 coordinates of the possibly
% matching points from pic_a and pic_b. Each row is a correspondence (e.g.
% row 42 of matches_a is a point that corresponds to row 42 of matches_b.

% 'Best_Fmatrix' is the 3x3 fundamental matrix
% 'inliers_a' and 'inliers_b' are the Mx2 corresponding points (some subset
% of 'matches_a' and 'matches_b') that are inliers with respect to
% Best_Fmatrix.

% For this section, use RANSAC to find the best fundamental matrix by
% randomly sample interest points. You would reuse
% estimate_fundamental_matrix() from part 2 of this assignment.

% If you are trying to produce an uncluttered visualization of epipolar
% lines, you may want to return no more than 30 points for either left or
% right images.

function [ Best_Fmatrix, inliers_a, inliers_b] = ransac_fundamental_matrix(matches_a, matches_b)


%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%

% Your ransac loop should contain a call to 'estimate_fundamental_matrix()'
% that you wrote for part II.

% set the varibles
[num, ~] = size(matches_a);
idx = 1 : 1 : num;
s = 8;
p = 0.99;
outlier_ratio = 0.5;
times = log(1 - p) / log(1-(1 - outlier_ratio) ^ 9);
inlier_persent = 0;
threshold = 0.05;
Best_Fmatrix = zeros(3, 3);
% do iterations on samples
for i = 1 : times
    x = randsample(idx, s);
    chosen_a = zeros(s, 2);
    chosen_b = zeros(s, 2);
    for j = 1 : s
        chosen_a(j, :) = matches_a(x(j), :);
        chosen_b(j, :) = matches_b(x(j), :);
    end
    matrix_f = estimate_fundamental_matrix(chosen_a, chosen_b);
    % check the accuracy of the model
    inlier_count = 0;
    for j = 1 : num
        distance = abs([matches_b(j, :), 1] * matrix_f * [matches_a(j, :), 1]');
        if distance <  threshold
            inlier_count = inlier_count + 1;
        end
    end
    if inlier_count / num > inlier_persent
        inlier_persent = inlier_count / num;
        Best_Fmatrix = matrix_f;
    end
end

% find the pairs
distances = zeros(num, 1);
for i = 1 : num
    distances(i) = abs([matches_b(i, :) 1] * Best_Fmatrix * [matches_a(i, :), 1]');
end
[B, oldIDX] = sort(distances);
len = min(length(find(B < threshold)), 30);
inliers_a = zeros(len, 2);
inliers_b = zeros(len, 2);
for i = 1 : len
    inliers_a(i, :) = matches_a(oldIDX(i), :);
    inliers_b(i, :) = matches_b(oldIDX(i), :);
end
end

