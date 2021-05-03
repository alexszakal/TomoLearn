clear all;
% The script generates a digital phantom with a rectangle in the center. The 
% attenuation of the rectangle is 400 HU everywhere.

P=uint16(zeros(1024,1024));
P(400:600, 200:800)=1400;
imshow(P)
caxis([0 max(max(P))])
title('rectangle.png')
figure(2);
plot(int32(P(512,:))-1000);
title('rectangle.png')
imwrite(P, 'rectangle.png', 'BitDepth', 16); % 16 bit

