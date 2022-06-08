function [x, y, mag, thresh, di]=SobelEdgeDetect(img, tlow, thigh)
[row, col, dim] = size(img);
Gx=[-1 0 1;-2 0 2;-1 0 1];
Gy=[1 2 1; 0 0 0; -1 -2 -1];
padsize=127;
padded=uint8(ones(row+2*padsize, col+2*padsize,1));
padded=padded.*127;
padded((1+padsize):(row+padsize),(1+padsize):(col+padsize),1)=img(:,:,1);
%W=double((1/((padsize*2+1)^2))*ones(padsize*2+1)); %kernel
outputimage = padded;
outputimage1 = padded;
outputimage2 = padded;
outputimage3 = padded;
threshold = padded;

for i=(1+padsize):(row+padsize) %loop over padded+ region
    for j=(1+padsize):(col+padsize)
        subregion=padded((i-1):(i+1),(j-1):(j+1));
        
        temp=double(subregion).*Gx; %kernel element-wise multiplication
        sumtemp=sum(temp(:));
        outputimage(i,j,1)=uint8(sumtemp);
        
        temp1=double(subregion).*Gy; %kernel element-wise multiplication
        sumtemp1=sum(temp1(:));
        outputimage1(i,j,1)=uint8(sumtemp1);

        outputimage2(i,j,1)=uint8(sqrt(sumtemp1^2+sumtemp^2)); %magnitude function
        
        %direction

        if rad2deg(atan2(sumtemp1,sumtemp))<0 %convert to two quadrants
            direc=360+rad2deg(atan2(sumtemp1,sumtemp));
            if direc>180
                direc=direc-180;
            end
        elseif rad2deg(atan2(sumtemp1,sumtemp))>180 
            direc=rad2deg(atan2(sumtemp1,sumtemp))-180;
        else
            direc=rad2deg(atan2(sumtemp1,sumtemp));
        end
        angles=[0 30 60 90 120 150 180];
        for l=1:6 %quantize
            if direc>=angles(l) && direc<=angles(l+1)
                direction=(angles(l+1)+angles(l))/2;
            end
        end
        outputimage3(i,j,1)=uint8(direction);
    end
end

for i=(1+padsize):(row+padsize) %loop over padded+ region
    for j=(1+padsize):(col+padsize)
        if outputimage2(i,j,1)<tlow %hysteresis
            threshold(i,j,1)=0;
            %outputimage3(i,j,1)=0;
        elseif outputimage2(i,j,1)<thigh %check for neighbors
            threshold(i,j,1)=0;
            %outputimage3(i,j,1)=0;
            for k=-1:1
                for p=-1:1
                    if outputimage2(i+k,j+p,1)>=thigh
                        threshold(i,j,1)=255;                        
                    end
                end
            end
        end
    end
end

              
x=outputimage((1+padsize):(row+padsize),(1+padsize):(col+padsize),1); %remove padding
y=outputimage1((1+padsize):(row+padsize),(1+padsize):(col+padsize),1);
mag=outputimage2((1+padsize):(row+padsize),(1+padsize):(col+padsize),1);
di=outputimage3((1+padsize):(row+padsize),(1+padsize):(col+padsize),1);
thresh=threshold((1+padsize):(row+padsize),(1+padsize):(col+padsize),1);
end