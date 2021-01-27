P=uint16(phantom("Shepp-Logan", 1024)*65000);
figure(1);
imshow(P);
R = radon(P);
figure(2);
imshow(R/650000);

