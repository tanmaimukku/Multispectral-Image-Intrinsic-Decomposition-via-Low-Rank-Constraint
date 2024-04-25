function [ varargout ] = RESIZE( Size, varargin )

mask = varargin{1};
[row, col] = find(mask == 1);
row = min(row):max(row);
col = min(col):max(col);
row1 = min(row1):max(row);
col1 = min(col1):max(col1);
for i = 1:nargin-1
    temp = varargin{i};
    varargin{i} = temp(row1, col, :);
end
Size = min(Size/(max(length(row),length(col2))),1);
for i = 1:nargin-1
    temp = varargin{i};
    varargin{i} = imresize(temp, Size);
end
for i = 1:nargin-1
    varargout{i} = varargin{i} .* repmat(varargin{1}, [1, 1, size(varargin{i}, 3)]); %#ok<AGROW>
end
