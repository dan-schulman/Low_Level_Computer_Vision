function [voting]=HoughTransformLineDetection(img)
[row, col, dim] = size(img);

for i=1:row %loop over padded+ region
    for j=1:col
        if img(i,j,1)>255
            img(i,j,1)=255;
        end
    end
end


dsize=100;
thetasize=30;
voting=zeros(dsize,thetasize);

maxd=sqrt(row^2+col^2);

d_array=-maxd:2*maxd/(dsize-1):maxd;
theta_array=-pi/2:pi/(thetasize-1):pi/2;

%store values for all pixels

for x=1:row
    for y=1:col
        if img(x,y,1)==0
            continue;
        end
        for theta=-pi/2:pi/(thetasize-1):pi/2
            d=x*cos(theta)+y*sin(theta);
            d_index=uint16(ceil((d+maxd)*50/maxd));
            theta_index=uint16(ceil((theta+pi/2)*30/pi));
            if d_index==0
                d_index=1;
            end
            if theta_index==0
                theta_index=1;
            end
            voting(d_index, theta_index)=voting(d_index, theta_index)+1;
        end
    end
end
figure;
pcolor(theta_array,d_array,voting);
shading flat;
title('Hough Transform');
xlabel('theta (radians)');
ylabel('d (pixels)');
colormap('hot');

 [r,c]=size(voting);
    
    
padsize=2; 
padded=(zeros(r+2*padsize, c+2*padsize));
sup=padded;
padded((1+padsize):(r+padsize),(1+padsize):(c+padsize))=voting;
for i=(1+padsize):(r+padsize) 
    for j=(1+padsize):(c+padsize)
        subregion=padded((i-2):(i+2),(j-2):(j+2));
        if max(subregion(:))==padded(i,j)
            %padded((i-2):(i+2),(j-2):(j+2))=zeros(5,5);
            %padded(i,j)=subregion(3,3);
            sup(i,j)=subregion(3,3);
        else
            padded(i,j)=0;
        end
    end
end
%voting2=padded((1+padsize):(r+padsize),(1+padsize):(c+padsize));
voting2=sup((1+padsize):(r+padsize),(1+padsize):(c+padsize));

%select top 3 candidates
counter=3;
maxx=[zeros(1,counter)];
index=zeros(2,counter);
while counter~=0
for i=1:r
    for j=1:c
        if voting2(i,j)>maxx(counter)
            maxx(counter)=voting(i,j);
            index(1,counter)=i;
            index(2,counter)=j;
            voting2(i,j)=0;
        end
    end   
end
counter=counter-1;
end
disp(maxx);
disp(index);
index(1,:)=(index(1,:)-1).*2*maxd/(dsize-1)-maxd;
index(2,:)=(index(2,:)-1).*pi/(thetasize-1)-pi/2;
disp(index);

figure;
imshow(img);

x=[0; index(1,1)*sin(index(2,1))];
y=[0; index(1,1)*cos(index(2,1))];

m = (diff(y)/diff(x));

minv = -1/m;


line([x(2) x(2)+600],[y(2) y(2)+600*minv],'Color','blue')
line([x(2) x(2)-600],[y(2) y(2)-600*minv],'Color','blue')
line(x,y)

x2=[0; index(1,2)*sin(index(2,2))];
y2=[0; index(1,2)*cos(index(2,2))];
m = (diff(y2)/diff(x2));
% Slope of new line
minv = -1/m;
line([x2(2) x2(2)+600],[y2(2) y2(2)+600*minv],'Color','blue')
line([x2(2) x2(2)-600],[y2(2) y2(2)-600*minv],'Color','blue')
line(x2,y2)

x3=[0; index(1,3)*sin(index(2,3))];
y3=[0; index(1,3)*cos(index(2,3))];

m = (diff(y3)/diff(x3));
% Slope of new line
minv = -1/m;
line([x3(2) x3(2)+600],[y3(2) y3(2)+600*minv],'Color','blue')
line([x3(2) x3(2)-600],[y3(2) y3(2)-600*minv],'Color','blue')
line(x3,y3)

end
