I = phantom('Modified Shepp-Logan',1024);
figure(1)
imshow(I,[]);
colorbar;


[R, xp] = radon(I);

% figure(2)
% imshow(R, []);
% colorbar;

[RI, H]=iradon(R,0:179,'linear','Ram-Lak', 0.7);

figure(3);
imshow(RI,[]);
colorbar;

figure(2);
plot(H);

figure(4);
plot(I(821,:));
hold on;
plot(RI(821,:));