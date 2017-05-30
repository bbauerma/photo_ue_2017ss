% =========================================================================
% PHOTOGRAMMETRIE VERTIEFUNG SS2017
% PROTOKOLL 1 - FUNDAMENTALMATRIX
%
% NAME - MNR
% ...
% Sebastian Flöry - 0826399
% =========================================================================

clc;
format short;

%=======================================================================
% TIEPOINTS IN BOTH IMAGES P1 AND P2
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

% =========================================================================
if length(p1) < 8
    error('At least 8 corresponding tiepoints are needed for the calculation of F')
end

if length(p1) ~= length(p2)
    error('Please provide the same amount of tiepoints in both images...')
end

camera_constant = 3400;    %px
image_width = 3000;        %px
image_height = 2000;       %px
origin = [0.5 -0.5];       %px

principal_point = [(image_width - 1) / 2., -(image_height - 1) / 2.];

% =========================================================================
% CALIBRATION MATRIX OF THE FIRST (C1) AND SECOND IMAGE (C2)
% =========================================================================
C1 = [
    1 0 -principal_point(1);
    0 1 -principal_point(2);
    0 0 -camera_constant
];

C2 = [
    1 0 -principal_point(1);
    0 1 -principal_point(2);
    0 0 -camera_constant
];

% =========================================================================
% CALCULATION OF THE FUNDAMENTALMATRIX
% =========================================================================
[m,n] = size(p1);

% initialize A which will store the kronecker product of corresponding
% tiepoints in each row; as we are using homogeneous coords the shape of
% A will be [nr_of_tiepoints, 3]
A = zeros(m,9);

for k = 1:m
    
    p1_hom_vec = [p1(k,2); p1(k,3); 1];
    p2_hom_vec = [p2(k,2); p2(k,3); 1];
    
    A(k,:) = transpose(kron(p2_hom_vec, p1_hom_vec));
    
end

AT = transpose(A);
AT_A = AT * A;

%SVD - SINGLE VALUE DECOMPOSITION
%fundamental matrix F is the eigenvector corresponding to the smallest
%eigenvalue in S; as the eigenvalues in S are sorted DESC the smallest
%is S[-1]; therefore the corresponding eigenvector is V(:,end)

[U,S,V] = svd(AT_A);
f = [V(:,end)];
F = reshape(f,[3,3]);

%as in general the determinatne of F won't be 0 we need to apply the SVD
%again; thistime we set the smalles eigenvalue in S to zero; by
%recalculating the fundamentalmatrix F0 using the adjusted eigenvalues we
%get the closest singular matrix to the original F

[U,S,V] = svd(F);

S(end,end) = 0;

F0 = U*S*V';

%caclulate the epipolar contradictions with both fundamentalmatrices
epi_con_F0 = 0;
epi_con_F = 0;

for k = 1:m
    
    p1_hom_vec = [p1(k,2); p1(k,3); 1];
    p2_hom_vec = [p2(k,2); p2(k,3); 1];
    
    epi_con_F0 = epi_con_F0 + transpose(p1_hom_vec)*F0*p2_hom_vec;
    epi_con_F = epi_con_F + transpose(p1_hom_vec)*F*p2_hom_vec;
end

avg_epi_con_F0 = epi_con_F0 / m;
avg_epi_con_F = epi_con_F / m;


% =========================================================================
% CALCULATION OF THE EPIPOLS (GAUSS EQUATIONS) - A*x = b
% F * e2 = 0
% F' * e1 = 0
% =========================================================================

%E1
A = [F0(1,1) F0(1,2); F0(2,1) F0(2,2); F0(3,1) F0(3,2)];
b = [-F0(1,3);-F0(2,3);-F0(3,3)];

%AT = transpose(A);
%ATA = AT * A;
%e1 = inv(ATA) * AT * b
e1 = A\b;
e1(3,1) = 1;

%E2 
F0T = transpose(F0);

A = [F0T(1,1) F0T(1,2);F0T(2,1) F0T(2,2);F0T(3,1) F0T(3,2)];
b = [-F0T(1,3);-F0T(2,3);-F0T(3,3)];

%AT = transpose(A);
%ATA = AT * A;
%e2 = inv(ATA) * AT * b;
e2 = A\b;
e2(3,1) = 1;


% =========================================================================
% RELATIVE ORIENTATION OF DEPENDENT IMAGES (FOLGEBILDANSCHLUSS)
% calculation of the exterior orientation based on the relative orientation
% of dependent images; therfore the exterior parameters of the first image
% are said to be zero; the parameters of the second image are relative to
% this one
% =========================================================================
R1 = eye(3);
Z1 = [0;0;0];

E = inv(transpose(C1))*F0*inv(C2);

[U,S,V] = svd(E);

Z2_1 = U(:,3);
Z2_2 = U(:,3) * (-1);

W = [0 -1 0; 1 0 0; 0 0 1];

R2_1 = U * W * V;
R2_2 = U * transpose(W) * V;

%there are four different solutions for the combination of R2 / Z2; we
%search for the one solution which indicates positive orientations of both
%images

if det(R2_1) < 0
    R2_1 = R2_1 * (-1);
end

if det(R2_2) < 0
    R2_2 = R2_2 * (-1);
end

p1_hom_vec = [p1(1,2); p1(1,3); 1];
p2_hom_vec = [p2(1,2); p2(1,3); 1];

%FIND THE CORRECT SOLUTION FOR R AND Z
for i = 1:4
    
    if i == 1
        Z2 = Z2_1;
        R2 = R2_1;
        
    elseif i == 2
        
        Z2 = Z2_1;
        R2 = R2_2;
        
    elseif i == 3
        
        Z2 = Z2_2;
        R2 = R2_1;
        
    elseif i == 4
        
        Z2 = Z2_2;
        R2 = R2_2;
   
    end
    
    S = C1 * p1_hom_vec;
    Q = R2 * C2 * p2_hom_vec;
    
    A = [S(1) -Q(1); S(2) -Q(2); S(3) -Q(3)];
    b = Z2;
    
    k = A\b;

    if all(k > 0)
        break
        
    end
end

% =========================================================================
% DISPLAY OUTPUT
% =========================================================================
disp('FUNDAMENTALMATRIX WITH DET!=0')
disp(F)
disp(' ')
disp('AVERAGE EPIPOLAR CONTRADIDCITON')
disp(avg_epi_con_F)
disp(' ')
disp('------------------------------------')
disp('FUNDAMENTALMATRIX WITH DET==0')
disp(F0)
disp(' ')
disp('AVERAGE EPIPOLAR CONTRADIDCITON')
disp(avg_epi_con_F0)
disp(' ')
disp(' ')
disp('============================================================')
disp('RELATIVE ORIENTATION OF DEPENDENT IMAGES (FOLGEBILDANSCHLUSS)')
disp('============================================================')
disp('[K1 ; K2]')
disp(k)
disp(' ')
disp('R1 - ROTATION MATRIX IMAGE 1 ')
disp(R1)
disp(' ')
disp('R2 - ROTATION MATRIX IMAGE 2 ')
disp(R2)
disp(' ')
disp('Z1 - PROJECTION CENTER OF FIRST IMAGE IN GLOBAL COORDINATES ')
disp(Z1)
disp(' ')
disp('Z2 - PROJECTION CENTER OF SECOND IMAGE IN GLOBAL COORDINATES ')
disp(Z2)
 