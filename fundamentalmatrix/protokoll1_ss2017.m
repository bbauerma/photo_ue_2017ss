

clc;
format short;
%=======================================================================
% test set used for controlling the results of the calculation of the
% fundamentalmatrix; points can be found in KRAUSS - PHOTOGRAMMETRIE BAND 1
% CHAPTER 4.3.1 - Page 219
% ======================================================================
p1=[
	54 1505.791 -934.224
	81 794.0793 -943.624
	85 1429.979 -1219.010
	129 831.493 -1154.156
	159 1722.855 -1221.836
	193 889.724 -1135.457
	215 1778.163 -1298.671
	237 1536.589 -901.604
	361 1049.607 -514.139
	407 1944.281 -1256.023
	468 1639.276 -1051.416
	510 1867.398 -716.761
];

p2=[
	54 1893.071 -718.891
	81 1217.832 -1371.529
	85 1708.045 -1159.236
	129 1159.775 -1541.060
	159 1960.698 -894.555
	193 1105.454 -1529.853
	215 1912.039 -981.543
	237 1594.205 -829.155
	361 1024.032 -1060.022
	407 1730.554 -915.923
	468 1227.926 -1121.613
	510 1505.277 -458.463
];

camera_constant = 3400;   %px

image_width = 3000;        %px
image_height = 2000;       %px
origin = [0.5 -0.5];     %px

principal_point = [(image_width - 1) / 2., -(image_height - 1) / 2.];

calibration_matrix = [
    1 0 -principal_point(1);
    0 1 -principal_point(2);
    0 0 -camera_constant;
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
epipolar_contradiction_F0 = 0;
epipolar_contradiction_F = 0;
for k = 1:m
    
    p1_hom_vec = [p1(k,2); p1(k,3); 1];
    p2_hom_vec = [p2(k,2); p2(k,3); 1];
    
    epipolar_contradiction_F0 = epipolar_contradiction_F0 + transpose(p1_hom_vec)*F0*p2_hom_vec;
    epipolar_contradiction_F = epipolar_contradiction_F + transpose(p1_hom_vec)*F*p2_hom_vec;
end

avg_epi_con_F0 = epipolar_contradiction_F0 / m 
avg_epi_con_F = epipolar_contradiction_F / m 
