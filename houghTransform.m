function [rho,theta,houghSpace] = houghTransform(theImage,thetaSampleFrequency)
 
    %Define the hough space
    theImage = flipud(theImage);
    [width,height] = size(theImage);
 
    rhoLimit = norm([width height]);
    rho = (-rhoLimit:1:rhoLimit);          
    theta = (0:thetaSampleFrequency:pi);
 
    numThetas = numel(theta);
    houghSpace = zeros(numel(rho),numThetas);
    %Find the "edge" pixels
    [xIndicies,yIndicies] = find(theImage);
 
    %Preallocate space for the accumulator array
    numEdgePixels = numel(xIndicies);
    accumulator = zeros(numEdgePixels,numThetas);
 
    cosine = (0:width-1)'*cos(theta); %Matrix Outerproduct  
    sine = (0:height-1)'*sin(theta); %Matrix Outerproduct
 
    accumulator((1:numEdgePixels),:) = cosine(xIndicies,:) + sine(yIndicies,:);
 
    %Scan over the thetas and bin the rhos 
    for i = (1:numThetas)
        houghSpace(:,i) = hist(accumulator(:,i),rho);
    end
    pcolor(theta,rho,houghSpace);
    shading flat;
    title('Hough Transform');
    xlabel('Theta (radians)');
    ylabel('Rho (pixels)');
    colormap('gray');
    
    [r,c]=size(houghSpace);
    
    
    padsize=2; 
padded=(zeros(r+2*padsize, c+2*padsize));
padded((1+padsize):(r+padsize),(1+padsize):(c+padsize))=houghSpace;
for i=(1+padsize):(r+padsize) 
    for j=(1+padsize):(c+padsize)
        subregion=padded((i-2):(i+2),(j-2):(j+2));
        if max(subregion(:))==padded(i,j)
            padded((i-2):(i+2),(j-2):(j+2))=zeros(5,5);
            padded(i,j)=subregion(3,3);
        else
            padded(i,j)=0;
        end
    end
end
houghSpace2=padded((1+padsize):(r+padsize),(1+padsize):(c+padsize));

%select top 3 candidates
counter=3;
maxx=[0 0 0];
index=zeros(2,3);
while counter~=0
for i=1:r
    for j=1:c
        if houghSpace2(i,j)>maxx(counter)
            maxx(counter)=houghSpace2(i,j);
            index(1,counter)=i;
            index(2,counter)=j;
            houghSpace2(i,j)=0;
        end
    end   
end
counter=counter-1;
end
disp(maxx);
disp(index);
index(1,:)=(index(1,:)-1).*2*rhoLimit/(r-1)-rhoLimit;
index(2,:)=(index(2,:)-1).*pi/(c-1);
disp(index);

figure;
theImage = flipud(theImage);
imshow(theImage);
y=[0; index(1,1)*sin(index(2,1))];
x=[0; index(1,1)*cos(index(2,1))];

m = (diff(y)/diff(x));
% Slope of new line
minv = -1/m;
line([x(2) x(2)+600],[y(2) y(2)+600*minv],'Color','blue')
line([x(2) x(2)-600],[y(2) y(2)-600*minv],'Color','blue')
line(x,y)

y2=[0; index(1,2)*sin(index(2,2))];
x2=[0; index(1,2)*cos(index(2,2))];
m = (diff(y2)/diff(x2));
% Slope of new line
minv = -1/m;
line([x2(2) x2(2)+600],[y2(2) y2(2)+600*minv],'Color','blue')
line([x2(2) x2(2)-600],[y2(2) y2(2)-600*minv],'Color','blue')
line(x2,y2)

y3=[0; index(1,3)*sin(index(2,3))];
x3=[0; index(1,3)*cos(index(2,3))];

m = (diff(y3)/diff(x3));
% Slope of new line
minv = -1/m;
line([x3(2) x3(2)+600],[y3(2) y3(2)+600*minv],'Color','blue')
line([x3(2) x3(2)-600],[y3(2) y3(2)-600*minv],'Color','blue')
line(x3,y3)
    
    
end