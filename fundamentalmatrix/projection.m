function [x] = projection(Z, R, C, X)
P = inv(C) * transpose(R) * [eye(3), -Z];
xh = P*[X', 1]';
x = xh/xh(3);
x = x(1:2);
end

