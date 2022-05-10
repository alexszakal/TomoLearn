clear all; close all;

figure(1);

%% PHANTOM
phantomSize=1024;
true_object = abs(phantom(phantomSize));
subplot(2,3,1); imshow(true_object, [0, max(true_object(:))]);
title("Original phantom"); 

%% MEASURED
angles = 1:1:180;
[meas_data, tVals] = radon(true_object, angles);
subplot(2,3,2); imshow(meas_data, [0, max(meas_data(:))]);
title("Measured"); 

rec = ones(size(true_object)); 

sinogram_ones = ones(size(meas_data));
sens = iradon(sinogram_ones,  angles, 'none', phantomSize);

for it = 1:110 
    
    forProject = radon(rec, angles);
    ratio = meas_data ./ ( forProject + 1e-5 );
    subplot(2,3,3); imshow(ratio', [0, max(ratio(:))]); title("Ratio");
    
    backProj_ratio = iradon(ratio,  angles, 'none', phantomSize);
    
    rec = rec .* backProj_ratio ./ sens;
    subplot(2,3,4); imshow(rec, [0, max(rec(:))]); title("Reconst. image");

    %cut along line 821
    subplot(2,3,5); plot(1:1:phantomSize, true_object(821,:), 1:1:phantomSize, rec(821,:) )
    title("iteration="+it);
    pause(0.5);
end
