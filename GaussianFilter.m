function final=GaussianFilter(img, sigma, n)
[row, col, dim] = size(img);
if dim>1
    red = double(img(:,:,1)); % Red channel
    green = double(img(:,:,2)); % Green channel
    blue = double(img(:,:,3)); % Blue channel
    gray = (red+green+blue)./3;
    img=gray; %convert to gray
end

%K=ones(n); %n=kernel size;

halfsize = floor(n/2);
[xx,yy] = meshgrid(-halfsize:halfsize, -halfsize:halfsize);
tmp = exp(-1/(2*sigma^2) * (xx.^2 + yy.^2));
tmp=tmp.*1./sum(tmp(:)); %normalize kernel
W=tmp;
figure;
surf(W);
colormap(jet);
disp(W);
title('Kernel Contour Plot');
%img = imread('racecar.tif'); 
padsize=floor(n/2); %padding to apply kernel
padded=uint8(zeros(row+2*padsize, col+2*padsize,1));
padded((1+padsize):(row+padsize),(1+padsize):(col+padsize),1)=img(:,:,1);
%W=double((1/((padsize*2+1)^2))*ones(padsize*2+1)); %kernel
outputimage = padded;
for i=(1+padsize):(row+padsize) %loop over padded+ region
    for j=(1+padsize):(col+padsize)
        subregion=padded((i-padsize):(i+padsize),(j-padsize):(j+padsize));
        temp=double(subregion).*W; %kernel element-wise multiplication
        sumtemp=sum(temp(:));
        outputimage(i,j,1)=uint8(sumtemp);
    end
end
final=outputimage((1+padsize):(row+padsize),(1+padsize):(col+padsize),1); %remove padding
end