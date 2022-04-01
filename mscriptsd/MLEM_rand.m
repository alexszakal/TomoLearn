clear all; close all;

figure(1);

%% PHANTOM
phantomSize=1024;
true_object = abs(phantom(phantomSize))/200;
subplot(2,3,1); imshow(true_object, [0, max(true_object(:))]);
title("Original phantom"); 

%% MEASURED
I0=1e3;
angles = 1:1:180;
[meas_data, tVals] = radon(true_object, angles);
meas_data=I0*exp(-1*meas_data);
meas_data=poissrnd(meas_data);
meas_data=-1*log(meas_data/I0);
%Filter the negative values from meas_data
meas_data = meas_data-(meas_data<0).*meas_data;
%figure(2);
%imshow(meas_data, [0, max(meas_data(:))]);

subplot(2,3,2); imshow(meas_data, [0, max(meas_data(:))]);
axis square;
title("Measured"); 

rec = ones(size(true_object)); 

sinogram_ones = ones(size(meas_data));
sens = iradon(sinogram_ones,  angles, 'none', phantomSize);

subplot(2,3,4); fbpRec=iradon(meas_data, angles, 'linear', 'Shepp-Logan', 0.6); imshow(fbpRec, [0, max(fbpRec(:))]);

for it = 1:200 
    
    forProject = radon(rec, angles);
    ratio = meas_data ./ ( forProject + 1e-5 );
    subplot(2,3,3); imshow(ratio', [0, max(ratio(:))]); title("Ratio");
    axis square;
    
    backProj_ratio = iradon(ratio,  angles, 'none', phantomSize);
    
    rec = rec .* backProj_ratio ./ sens;
    subplot(2,3,5); imshow(rec, [0, max(rec(:))]); title("Reconst. image");

    %cut along line 821
    subplot(2,3,6); plot(1:1:phantomSize, true_object(821,:), 1:1:phantomSize, rec(821,:) )
    title("iteration="+it);
    pause(0.5);
end

 
pause(10)
