

clc;
format short;
%=======================================================================
% test set used for controlling the results of the calculation of the
% fundamentalmatrix; points can be found in KRAUSS - PHOTOGRAMMETRIE BAND 1
% CHAPTER 4.3.1 - Page 219
% ======================================================================
p1=[
    1 93.176 5.890;
    2 -27.403 6.672;
    3 83.951 107.422;
    4 -11.659 101.544;
    5 110.326 -97.800;
    6 -12.653 -87.645;
    7 37.872 40.969;
    8 41.503 -37.085
];

p2=[
    1 6.072 5.176;
    2 -112.842 1.121;
    3 -4.872 105.029;
    4 -99.298 95.206;
    5 34.333 -99.522;
    6 -96.127 -93.761;
    7 -48.306 37.862;
    8 -42.191 -40.138
];

camera_constant = 152.67;   %mm

image_width = 3000;        %px
image_height = 2000;       %px
origin = [0.5 -0.5];     %px

principal_point = [(img_width - 1) / 2., -(img_height - 1) / 2.];

calibration_matrix = [
    1 0 -principal_point(1);
    0 1 -principal_point(2);
    0 0 -camera_const;
];

[m,n] = size(p1);

% A contains the kronecker product of P2 and P1 in each row
A = zeros(m,9);

for k = 1:m
    
    p1_hom_vec = [p1(k,2); p1(k,3); 1];
    p2_hom_vec = [p2(k,2); p2(k,3); 1];
    
    A(k,:) = transpose(kron(p2_hom_vec, p1_hom_vec));
    
end

AT = transpose(A);
AT_A = AT * A;

%SVD - SINGLE VALUE DECOMPOSITION
%fundamental matrix f is the eigenvector corresponding to the smallest
%eigenvalue in s; as the eigenvalues in s are sorted DESC the smallest
%is s[-1]; the resultung fundamentalmatrix f in general won't have
%det(F) == 0

[U,S,V] = svd(AT_A);

f = [V(:,end)];

F = reshape(f,[3,3]);
disp('FUNDAMENTALMATRIX DET != 0')
F = F / F(end,end)

[U,S,V] = svd(F);

S(end,end) = 0;
 
disp('FUNDAMENTALMATRIX DET = 0')
F0 = U*S*V'

%check if F / F1 are correct:
% for k = 1:m
%     
%     p1_hom_vec = [p1(k,2); p1(k,3); 1];
%     p2_hom_vec = [p2(k,2); p2(k,3); 1];
%     
%     epipolar_contradiction = transpose(p1_hom_vec)*F0*p2_hom_vec
%     epipolar_contradiction = transpose(p1_hom_vec)*F*p2_hom_vec
% end
