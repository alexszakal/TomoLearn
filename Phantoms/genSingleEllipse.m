x0=0.2;
y0=0.3;
a=0.3;
b=0.45;
A=1.0;
phi=30;

numPoints=1024;
I = phantom([A,a,b,x0,y0,phi],numPoints);
imwrite(I, 'singleEllipse.png', 'bitDepth', 16);
figure(1)
imshow(I);

theta=90;
[R, xp] = radon(I, theta);
aParam = sqrt( (a/2*numPoints)^2 * cosd(theta-phi)^2 + (b/2*numPoints)^2 * sind(theta-phi)^2 );
s=sqrt( (x0/2*numPoints)^2 + (y0/2*numPoints)^2 );
if x0==0 && y0==0
    gamma=0;
else
    gamma=atan2d(y0,x0);
end

aR=R*0;
for i=1:size(xp)
    if abs(xp(i)-s*cosd(gamma-theta))  < aParam
        aR(i) = 2*A*(a*numPoints/2)*(b*numPoints/2)/aParam^2* sqrt(aParam^2-(xp(i)-s*cosd(gamma-theta))^2);
    end
end


figure(2)
plot(xp, R, xp, aR*1.00);





