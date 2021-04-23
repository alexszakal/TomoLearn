clear all;

[P,E]=phantom('Modified Shepp-Logan', 1024);
E(:,1)=[2., -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
maxInHU=1000; %Attenuation of bone in HU
P=uint16(phantom(E, 1024)*maxInHU/2*2);
imshow(P);
caxis([0 maxInHU])
title('Shepp-Logan.png')
figure(2);
plot(P(822,:));
title('Shepp-Logan.png')
imwrite(P, 'SheppLogan.png', 'BitDepth', 16); % 16 bit

P(1000:1024,:)=[];
P(1:24,:)=[];
P(:,1:130)=[];
imwrite(P, 'SheppLogan_asymmetric.png', 'BitDepth', 16); %16 bit


P2=uint16(phantom('Modified Shepp-Logan', 1024)*maxInHU*2);
figure(3);
imshow(P2);
ax=gca;
caxis(ax, [0 maxInHU])
title('ModifiedSheppLogan.png')
figure(4);
plot(P2(822,:));
title('ModifiedSheppLogan.png')
imwrite(P2, 'ModifiedSheppLogan.png', 'BitDepth', 16);  %16/15 bit

P2(1000:1024,:)=[];
P2(1:24,:)=[];
P2(:,1:130)=[];
%figure(5);
%imshow(P2);
imwrite(P2, 'ModifiedSheppLogan_asymmetric.png', 'BitDepth', 16); % 16/15 bit
