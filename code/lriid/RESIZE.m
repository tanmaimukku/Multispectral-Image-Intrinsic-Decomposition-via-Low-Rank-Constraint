function [ varargout ] = RESIZE( Size, varargin )

mask = varargin{1};

[block1, block2] = find(mask == 1);
block1 = min(block1):max(block1);
block2 = min(block2):max(block2);
[block11, block12] = find(mask == 1);
block11 = min(block11):max(block11);
block12 = min(block12):max(block12);
for i = 1:nargin-1
    temp = varargin{i};
    varargin{i} = temp(block11, block2, :);
end
for i = 1:nargin-1
    temp = varargin{i};
    varargin{i} = imresize(temp, Size);
end
for i = 1:nargin-1
    varargout{i}=varargin{i}.*repmat(varargin{1},[1, 1, size(varargin{i},3)]);
end
