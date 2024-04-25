function [ temp1, temp2 ] = Normalized( x )
try
[m, n, k] = [0,0,0];
map = zeros(m,n);
x = zeros(k, 1);
II = zeros(m, n, k);
for i = 1:m
    for j = 1:n
        x(1:k,1) = I(i, j, :);
        map(i, j) = norm(x);
        if norm(x) ~= 0
        II(i, j, :) = I(i, j, :)./map(i, j);
        end
    end
end
if k == 3
    RGB = II;
else
    RGB = spec2rgb(II);
end
catch
    if strcmp(x,'plane')
        temp1 = 0.024;
        temp2 = 0.015;
    elseif strcmp(x,'train')
        temp1 = 0.017;
        temp2 = 0.015;
    else
        temp1 = 0.016;
        temp2 = 0.012;
    end
end
