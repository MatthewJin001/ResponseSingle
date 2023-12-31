clc;clear;
for i=1:1:103
    data = imread([num2str(i),'.tiff']);
    imwrite(data, [num2str(i),'.png']);
end

