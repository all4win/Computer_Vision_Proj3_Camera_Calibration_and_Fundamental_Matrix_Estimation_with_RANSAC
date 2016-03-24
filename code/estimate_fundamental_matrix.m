% Fundamental Matrix Stencil Code
% CS 4495 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Returns the camera center matrix for a given projection matrix

% 'Points_a' is nx2 matrix of 2D coordinate of points on Image A
% 'Points_b' is nx2 matrix of 2D coordinate of points on Image B
% 'F_matrix' is 3x3 fundamental matrix

% Try to implement this function as efficiently as possible. It will be
% called repeatly for part III of the project

function [ F_matrix ] = estimate_fundamental_matrix(Points_a,Points_b)

%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%

% normalization matrix
S_a = diag([1, 1, 1]);
S_b = diag([1, 1, 1]);
C_a = diag([1, 1, 1]);
C_b = diag([1, 1, 1]);
ca = mean(Points_a);
cb = mean(Points_b);
C_a(7: 8) = -ca(1 : 2);
C_b(7: 8) = -cb(1 : 2);
std_cof = 2 ^ (0.5);
std_a = std(Points_a(:,1) - ca(1)) + std(Points_a(:,2) - ca(2));
std_b = std(Points_b(:,1) - cb(1)) + std(Points_b(:,2) - cb(2));
sa = std_cof / std_a;
sb = std_cof / std_b;
S_a(1) = sa;
S_a(5) = sa;
S_b(1) = sb;
S_b(5) = sb;
T_a = S_a * C_a;
T_b = S_b * C_b;
Points_a(:, 3) = 1;
Points_b(:, 3) = 1;
% normalize the points
Points_a = (T_a * Points_a')';
Points_b = (T_b * Points_b')';


row = size(Points_a, 1);
A = zeros(row, 8);
B = -ones(row, 1);
for i = 1 : row
    A(i, :) = [Points_a(i, 1) * Points_b(i, 1), Points_a(i, 2) * Points_b(i, 1), Points_b(i, 1), Points_a(i, 1) * Points_b(i, 2), Points_a(i, 2) * Points_b(i, 2), Points_b(i, 2), Points_a(i, 1), Points_a(i, 2)];
end
F = A \ B;
F(9) = 1;
F_matrix = reshape(F, [3, 3])';
[U, S, V] = svd(F_matrix);
S(3, 1:3) = zeros(1, 3);
F_matrix = U * S * V';

% normalize the fundamental matrix
F_matrix = T_b' * F_matrix * T_a;
end

