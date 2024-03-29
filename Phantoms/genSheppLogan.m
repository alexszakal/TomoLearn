clear all;
% Script to generate the digital phantoms for reconstruction
% The phantoms created by the standard Matlab functions then divided by 2.5
% to get realistic values in [1/cm] units. Then the mu values are converted to
% Hounsfield units (mu_air ~= 0 was used). 
% Since the png can handle only positive values and the vacuum has -1000 HU value,
% 1000 was added to the values. 

[P,E]=phantom('Shepp-Logan', 1024);
E(:,1)=[2., -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
mu_water=0.2; %[1/cm]
%P=uint16(1000*(phantom(E, 1024)/5-mu_water)/mu_water + 1000);
P=uint16(1000*(P/2.5-mu_water)/mu_water + 1000);
imshow(P);
caxis([0 max(max(P))])
title('Shepp-Logan.png')
figure(2);
plot(int32(P(822,:))-1000);
title('Shepp-Logan.png')
imwrite(P, 'SheppLogan_HU.png', 'BitDepth', 16); % 16 bit

P(1000:1024,:)=[];
P(1:24,:)=[];
P(:,1:130)=[];
imwrite(P, 'SheppLogan_asymmetric_HU.png', 'BitDepth', 16); %16 bit


P2=uint16(1000*(phantom('Modified Shepp-Logan', 1024)/2.5-mu_water)/mu_water + 1000);
figure(3);
imshow(P2);
ax=gca;
caxis([0 max(max(P2))])
title('ModifiedSheppLogan.png')
figure(4);
plot(int32(P2(822,:))-1000);
title('ModifiedSheppLogan.png')
imwrite(P2, 'ModifiedSheppLogan_HU.png', 'BitDepth', 16);  %16/15 bit

P2(1000:1024,:)=[];
P2(1:24,:)=[];
P2(:,1:130)=[];
%figure(5);
%imshow(P2);
imwrite(P2, 'ModifiedSheppLogan_asymmetric_HU.png', 'BitDepth', 16); % 16/15 bit
