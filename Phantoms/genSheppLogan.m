[P,E]=phantom('Modified Shepp-Logan', 1024);
E(:,1)=[2., -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
P=uint16(phantom(E, 1024)*2^15);
imshow(P);
figure(2);
plot(P(822,:));
imwrite(P, 'SheppLogan.png', 'BitDepth', 16);

P(1000:1024,:)=[];
P(1:24,:)=[];
P(:,1:130)=[];
%figure(5);
%imshow(P);
imwrite(P, 'SheppLogan_asymmetric.png', 'BitDepth', 16);


P2=uint16(phantom('Modified Shepp-Logan', 1024)*2^16);
figure(3);
imshow(P2);
figure(4);
plot(P2(822,:));
imwrite(P2, 'ModifiedSheppLogan.png', 'BitDepth', 16);

P2(1000:1024,:)=[];
P2(1:24,:)=[];
P2(:,1:130)=[];
%figure(5);
%imshow(P2);
imwrite(P2, 'ModifiedSheppLogan_asymmetric.png', 'BitDepth', 16);
