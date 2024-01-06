clear, clc, close all
%% 1
% Load the hyperspectral images, crop them into a subregion of 900x1000px and save as new files
% hcube = hypercube('EL_GEO_120243.img','EL_GEO_120243.hdr');
% 
% row = 1550:2450;
% column = 4600:5600;
% newhcube = cropData(hcube,row,column,':');
% enviwrite(newhcube,'EL_GEO_120243_cropped');
% 
% hcube = hypercube('EL_GEO_082915.img','EL_GEO_082915.hdr');
% 
% row = 4600:5500;
% column = 4400:5400;
% newhcube = cropData(hcube,row,column,':');
% enviwrite(newhcube,'EL_GEO_082915_cropped');
% 
% row = 3700:4600;
% column = 3650:4650;
% newhcube = cropData(hcube,row,column,':');
% enviwrite(newhcube,'EL_GEO_082915_cropped_');

% Load the cropped hyperspectral images
fasano = hypercube('EL_GEO_120243_cropped.dat','EL_GEO_120243_cropped.hdr');
grottaglie = hypercube('EL_GEO_082915_cropped.dat','EL_GEO_082915_cropped.hdr');
%% 2
% Data analysis
sF=sign(fasano.DataCube);
minF=min(fasano.DataCube,[],'all')
maxF=max(fasano.DataCube,[],'all')
ipositifF=sum(sF(:)==1)
inegatifF=sum(sF(:)==-1)

sG=sign(grottaglie.DataCube);
minG=min(grottaglie.DataCube,[],'all')
maxG=max(grottaglie.DataCube,[],'all')
ipositifG=sum(sG(:)==1)
inegatifG=sum(sG(:)==-1)

% Clipping the values in the range [0, 1]
fasanoR = max(min(fasano.DataCube, 1), 0);
grottaglieR = max(min(grottaglie.DataCube, 1), 0);

fasano = assignData(fasano,':',':',':',fasanoR);
grottaglie = assignData(grottaglie,':',':',':',grottaglieR);

% Preprocess data (removing the last band (band 48) since it has all zero values in the 'grottaglie' hypercube)
bandsToRemove = [];
nBands = fasano.Metadata.Bands;
for i=1:nBands
    if nnz(grottaglie.DataCube(:,:,i))==0 || nnz(fasano.DataCube(:,:,i))==0
        bandsToRemove = [bandsToRemove, i];
    end
end
bandsToRemove
for i=1:length(bandsToRemove)
    bandToRemove = bandsToRemove(i);
    fasano = removeBands(fasano,'BandNumber',bandToRemove);
    grottaglie = removeBands(grottaglie,'BandNumber',bandToRemove);
end

% Estimate RGB color images from the data cubes
rgbImg_fasano = colorize(fasano,'Method','RGB','ContrastStretching',true);
figure
imagesc(rgbImg_fasano)
title('RGB Image of Fasano Cropped Area')
axis off

rgbImg_grottaglie = colorize(grottaglie,'Method','RGB','ContrastStretching',true);
figure
imagesc(rgbImg_grottaglie)
title('RGB Image of Grottaglie Cropped Area (Before Normalization)')
axis off
%% 3
% Use the imageSegmenter tool to create a mask for selecting some area corresponding to olive crowns in 'fasano' image 
% imageSegmenter(fasano.DataCube(:,:,47))
maskF = segmentMaskFasano(fasano.DataCube(:,:,47));
structBoundaries = bwboundaries(maskF);
figure, imagesc(rgbImg_fasano)
title('Selected Crowns (Fasano)')
axis off
hold on
for i=1:length(structBoundaries)
   xy=structBoundaries{i};
   x = xy(:, 2);
   y = xy(:, 1);
   plot(x, y, 'LineWidth', 2);
   drawnow 
end
hold off

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area 
[row, col] = find(maskF); % row and col where maskF is different from 0
reflectanceOliveM = [];
for i=1:length(row)
    reflectance = fasano.DataCube(row(i),col(i),:);
    reflectanceOliveM = [reflectanceOliveM, reshape(reflectance, [fasano.Metadata.Bands,1])];
end
reflectanceOlive = mean(reflectanceOliveM, 2);

% Use the imageSegmenter tool to create a mask for selecting some area corresponding to olive crowns in 'grottaglie' images
% imageSegmenter(grottaglie.DataCube(:,:,47))
maskG = segmentMaskGrottaglie(grottaglie.DataCube(:,:,47));
structBoundaries = bwboundaries(maskG);
figure, imagesc(rgbImg_grottaglie)
title('Selected Crowns (Grottaglie)')
axis off
hold on
for i=1:length(structBoundaries)
   xy=structBoundaries{i};
   x = xy(:, 2);
   y = xy(:, 1);
   plot(x, y, 'LineWidth', 2);
   drawnow 
end
hold off

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area 
[row, col] = find(maskG);
reflectanceOliveNotNormalizedM = [];
for i=1:length(row)
    reflectance = grottaglie.DataCube(row(i),col(i),:);
    reflectanceOliveNotNormalizedM = [reflectanceOliveNotNormalizedM, reshape(reflectance, [grottaglie.Metadata.Bands,1])];
end
reflectanceOliveNotNormalized = mean(reflectanceOliveNotNormalizedM, 2);

figure
hold on
plot(fasano.Wavelength(:),reflectanceOlive,'LineWidth',2)
plot(grottaglie.Wavelength(:),reflectanceOliveNotNormalized,'LineWidth',2)
axis tight
box on
title('Spectral Signature of Olive Tree Crown')
xlabel('Wavelength (\mum)')
ylabel('Reflectance')
legend('Actual','Not Normalized','Location','best')
title(legend,'Image')
hold off
%% 4
% Normalize the miscalibrated images using Histogram Matching
grottaglieN_HM = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);

for b=1:grottaglie.Metadata.Bands
    grottaglieN_HM(:,:,b) = imhistmatch(grottaglie.DataCube(:,:,b), fasano.DataCube(:,:,b)); 
end

% Histogram Matching using lut
% grottaglieN = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);
% for b=1:grottaglie.Metadata.Bands
%     h=imhist(grottaglie.DataCube(:,:,b))./numel(grottaglie.DataCube(:,:,b));
%     cdf_g=cumsum(h);
%     r=imhist(fasano.DataCube(:,:,b))./numel(fasano.DataCube(:,:,b));
%     cdf_f=cumsum(r);
%     LUT=zeros(256,1);
%     gj=0;
%     for g=0:255
%         while gj<255 && cdf_f(gj+1) < cdf_g(g+1)
%             gj=gj+1;
%         end
%         LUT(g+1)=gj;
%     end
%     LUT=uint8(LUT);
%     grottaglieN(:,:,b) = single(intlut(uint8(255 * grottaglie.DataCube(:,:,b)),LUT)) / 255;
% end

grottaglieN_HM = assignData(grottaglie,':',':',':',grottaglieN_HM);
 
rgbImg_grottaglieN_HM = colorize(grottaglieN_HM,'Method','RGB',"ContrastStretching",true);
figure
imagesc(rgbImg_grottaglieN_HM)
axis off
title('RGB Image of Grottaglie After Normalization (Histogram Matching)')

% Clipping the values in the range [0, 1]
grottaglieNR_HM = max(min(grottaglieN_HM.DataCube, 1), 0);
grottaglieN_HM = assignData(grottaglieN_HM,':',':',':',grottaglieNR_HM);

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area after normalization
[row, col] = find(maskG);
reflectanceOliveNormalizedHMM = [];
for i=1:length(row)
    reflectance = grottaglieN_HM.DataCube(row(i),col(i),:);
    reflectanceOliveNormalizedHMM = [reflectanceOliveNormalizedHMM, reshape(reflectance, [grottaglieN_HM.Metadata.Bands,1])];
end
reflectanceOliveNormalizedHM = mean(reflectanceOliveNormalizedHMM, 2);
%% 5
% Normalize the miscalibrated images using Minimum-Maximum normalization
grottaglieN_MM = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);

for band=1:grottaglie.Metadata.Bands
    min_f = min(fasano.DataCube(:,:,band),[],'all');
    max_f = max(fasano.DataCube(:,:,band),[],'all');
    min_g = min(grottaglie.DataCube(:,:,band),[],'all');
    max_g = max(grottaglie.DataCube(:,:,band),[],'all');
    
    a = (max_f - min_f) / (max_g - min_g);
    b = min_f - a * min_g;
    
    grottaglieN_MM(:,:,band) = grottaglie.DataCube(:,:,band) * a + b;
end

grottaglieN_MM = assignData(grottaglie,':',':',':',grottaglieN_MM);
 
rgbImg_grottaglieN_MM = colorize(grottaglieN_MM,'Method','RGB',"ContrastStretching",true);
figure
imagesc(rgbImg_grottaglieN_MM)
axis off
title('RGB Image of Grottaglie After Normalization (Minimum-Maximum)')

% Clipping the values in the range [0, 1]
grottaglieNR_MM = max(min(grottaglieN_MM.DataCube, 1), 0);
grottaglieN_MM = assignData(grottaglieN_MM,':',':',':',grottaglieNR_MM);

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area after normalization
[row, col] = find(maskG);
reflectanceOliveNormalizedMMM = [];
for i=1:length(row)
    reflectance = grottaglieN_MM.DataCube(row(i),col(i),:);
    reflectanceOliveNormalizedMMM = [reflectanceOliveNormalizedMMM, reshape(reflectance, [grottaglieN_MM.Metadata.Bands,1])];
end
reflectanceOliveNormalizedMM = mean(reflectanceOliveNormalizedMMM, 2);
%% 6
% Normalize the miscalibrated images using Mean-Standard Deviation normalization
grottaglieN_MS = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);

for band=1:grottaglie.Metadata.Bands
    mean_f = mean(fasano.DataCube(:,:,band), 'all');
    std_f = std(fasano.DataCube(:,:,band), 0, 'all');
    mean_g = mean(grottaglie.DataCube(:,:,band), 'all');
    std_g = std(grottaglie.DataCube(:,:,band), 0, 'all');
    
    a = std_f / std_g;
    b = mean_f - a * mean_g;
    
    grottaglieN_MS(:,:,band) = grottaglie.DataCube(:,:,band) * a + b;
end

grottaglieN_MS = assignData(grottaglie,':',':',':',grottaglieN_MS);
 
rgbImg_grottaglieN_MS = colorize(grottaglieN_MS,'Method','RGB',"ContrastStretching",true);
figure
imagesc(rgbImg_grottaglieN_MS)
axis off
title('RGB Image of Grottaglie After Normalization (Mean-Standard Deviation)')

% Clipping the values in the range [0, 1]
grottaglieNR_MS = max(min(grottaglieN_MS.DataCube, 1), 0);
grottaglieN_MS = assignData(grottaglieN_MS,':',':',':',grottaglieNR_MS);

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area after normalization
[row, col] = find(maskG);
reflectanceOliveNormalizedMSM = [];
for i=1:length(row)
    reflectance = grottaglieN_MS.DataCube(row(i),col(i),:);
    reflectanceOliveNormalizedMSM = [reflectanceOliveNormalizedMSM, reshape(reflectance, [grottaglieN_MS.Metadata.Bands,1])];
end
reflectanceOliveNormalizedMS = mean(reflectanceOliveNormalizedMSM, 2);
%% 7
% Normalize the miscalibrated images using Simple Regression normalization
grottaglieN_SR = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);

for band=1:grottaglie.Metadata.Bands
    mean_f = mean(fasano.DataCube(:,:,band), 'all');
    mean_g = mean(grottaglie.DataCube(:,:,band), 'all');
    var_g = var(grottaglie.DataCube(:,:,band), 0, 'all');
    cov_g_f = cov(grottaglie.DataCube(:,:,band), fasano.DataCube(:,:,band));
     
    a = cov_g_f(1, 2) / var_g;
    b = mean_f - a * mean_g;
    
    grottaglieN_SR(:,:,band) = grottaglie.DataCube(:,:,band) * a + b;
end

grottaglieN_SR = assignData(grottaglie,':',':',':',grottaglieN_SR);
 
rgbImg_grottaglieN_SR = colorize(grottaglieN_SR,'Method','RGB',"ContrastStretching",true);
figure
imagesc(rgbImg_grottaglieN_SR)
axis off
title('RGB Image of Grottaglie After Normalization (Simple Regression)')

% Clipping the values in the range [0, 1]
grottaglieNR_SR = max(min(grottaglieN_SR.DataCube, 1), 0);
grottaglieN_SR = assignData(grottaglieN_SR,':',':',':',grottaglieNR_SR);

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area after normalization
[row, col] = find(maskG);
reflectanceOliveNormalizedSRM = [];
for i=1:length(row)
    reflectance = grottaglieN_SR.DataCube(row(i),col(i),:);
    reflectanceOliveNormalizedSRM = [reflectanceOliveNormalizedSRM, reshape(reflectance, [grottaglieN_SR.Metadata.Bands,1])];
end
reflectanceOliveNormalizedSR = mean(reflectanceOliveNormalizedSRM, 2);
%% 8
% Normalize the miscalibrated images using Pseudo Invariant Features normalization

% Load, preprocess and use the other cropped area of Grottaglie image data for computing the normalization coefficients
grottaglie_ = hypercube('EL_GEO_082915_cropped_.dat','EL_GEO_082915_cropped_.hdr');

grottaglieR_= max(min(grottaglie_.DataCube, 1), 0);
grottaglie_ = assignData(grottaglie_,':',':',':',grottaglieR_);

for i=1:length(bandsToRemove)
    bandToRemove = bandsToRemove(i);
    grottaglie_ = removeBands(grottaglie_,'BandNumber',bandToRemove);
end

rgbImg_grottaglie_ = colorize(grottaglie_,'Method','RGB','ContrastStretching',true);

maskG_ = segmentMaskGrottaglie_(grottaglie_.DataCube(:,:,47));
structBoundaries = bwboundaries(maskG_);
figure, imagesc(rgbImg_grottaglie_)
title('Selected PIFs for Image Normalization')
axis off
hold on
for i=1:length(structBoundaries)
   xy=structBoundaries{i};
   x = xy(:, 2);
   y = xy(:, 1);
   plot(x, y, 'LineWidth', 2);
   drawnow 
end
hold off

[row, col] = find(maskG_);
reflectanceOliveNotNormalizedM_ = [];
for i=1:length(row)
    reflectance = grottaglie_.DataCube(row(i),col(i),:);
    reflectanceOliveNotNormalizedM_ = [reflectanceOliveNotNormalizedM_, reshape(reflectance, [grottaglie_.Metadata.Bands,1])];
end
reflectanceOliveNotNormalized_ = mean(reflectanceOliveNotNormalizedM_, 2);

grottaglieN_PIF = zeros(grottaglie.Metadata.Height, grottaglie.Metadata.Width, grottaglie.Metadata.Bands);

mean_f = reflectanceOlive;
std_f = std(reflectanceOliveM, 0, 2);

mean_g = reflectanceOliveNotNormalized_;
std_g = std(reflectanceOliveNotNormalizedM_, 0, 2);

for band=1:grottaglie.Metadata.Bands
    a = std_f(band) / std_g(band);
    b = mean_f(band) - a * mean_g(band);
    
    % Normalization of the image using the normalization coeffients computed considering the other one
    grottaglieN_PIF(:,:,band) = grottaglie.DataCube(:,:,band) * a + b;
end

grottaglieN_PIF = assignData(grottaglie,':',':',':',grottaglieN_PIF);
 
rgbImg_grottaglieN_PIF = colorize(grottaglieN_PIF,'Method','RGB',"ContrastStretching",true);
figure
imagesc(rgbImg_grottaglieN_PIF)
axis off
title('RGB Image of Grottaglie After Normalization (Pseudo Invariant Features)')

% Clipping the values in the range [0, 1]
grottaglieNR_PIF = max(min(grottaglieN_PIF.DataCube, 1), 0);
grottaglieN_PIF = assignData(grottaglieN_PIF,':',':',':',grottaglieNR_PIF);

% Estimate the spectral signature of olive tree crown computing the mean of the pixels values for each band inside the masked area after normalization
[row, col] = find(maskG);
reflectanceOliveNormalizedPIFM = [];
for i=1:length(row)
    reflectance = grottaglieN_PIF.DataCube(row(i),col(i),:);
    reflectanceOliveNormalizedPIFM = [reflectanceOliveNormalizedPIFM, reshape(reflectance, [grottaglieN_PIF.Metadata.Bands,1])];
end
reflectanceOliveNormalizedPIF = mean(reflectanceOliveNormalizedPIFM, 2);
%% 9
% Plot the spectral signatures to test the normalization
figure
hold on
plot(fasano.Wavelength(:),reflectanceOlive,'LineWidth',2)
plot(grottaglieN_HM.Wavelength(:),reflectanceOliveNormalizedHM,'LineWidth',2)
plot(grottaglieN_MM.Wavelength(:),reflectanceOliveNormalizedMM,'LineWidth',2)
plot(grottaglieN_MS.Wavelength(:),reflectanceOliveNormalizedMS,'LineWidth',2)
plot(grottaglieN_SR.Wavelength(:),reflectanceOliveNormalizedSR,'LineWidth',2)
plot(grottaglieN_PIF.Wavelength(:),reflectanceOliveNormalizedPIF,'LineWidth',2)
axis tight
box on
title('Spectral Signature of Olive Tree Crown')
xlabel('Wavelength (\mum)')
ylabel('Reflectance')
legend('Calibrated','Normalized (Histogram Matching)','Normalized (Minimum-Maximum)','Normalized (Mean-Standard Deviation)','Normalized (Simple Regression)','Normalized (Pseudo Invariant Features)','Location','best')
title(legend,'Image')
hold off

% Compute the absolute reflectance differences for each bands between the actual olive tree crown spectral signature and the ones obtained with the different normalization methods
ASSD_HM = [];
ASSD_MM = [];
ASSD_MS = [];
ASSD_SR = [];
ASSD_PIF = [];

for band=1:fasano.Metadata.Bands
    ASSD_HM = [ASSD_HM, abs(reflectanceOlive(band)-reflectanceOliveNormalizedHM(band))]; 
    ASSD_MM = [ASSD_MM, abs(reflectanceOlive(band)-reflectanceOliveNormalizedMM(band))];
    ASSD_MS = [ASSD_MS, abs(reflectanceOlive(band)-reflectanceOliveNormalizedMS(band))];
    ASSD_SR = [ASSD_SR, abs(reflectanceOlive(band)-reflectanceOliveNormalizedSR(band))];
    ASSD_PIF = [ASSD_PIF, abs(reflectanceOlive(band)-reflectanceOliveNormalizedPIF(band))];
end

figure
hold on
plot(grottaglieN_HM.Wavelength(:),ASSD_HM,'LineWidth',2)
plot(grottaglieN_MM.Wavelength(:),ASSD_MM,'LineWidth',2)
plot(grottaglieN_MS.Wavelength(:),ASSD_MS,'LineWidth',2)
plot(grottaglieN_SR.Wavelength(:),ASSD_SR,'LineWidth',2)
plot(grottaglieN_PIF.Wavelength(:),ASSD_PIF,'LineWidth',2)
axis tight
box on
title('Absolute Spectral Signature Differences')
xlabel('Wavelength (\mum)')
ylabel('Absolute Reflectance Difference')
legend('Histogram Matching','Minimum-Maximum','Mean-Standard Deviation','Simple Regression','Pseudo Invariant Features','Location','best')
title(legend,'Image')
hold off

% Compute the mean absolute spectral signuture differences to test the normalization
MASSD_HM = mean(ASSD_HM);
MASSD_MM = mean(ASSD_MM);
MASSD_MS = mean(ASSD_MS);
MASSD_SR = mean(ASSD_SR);
MASSD_PIF = mean(ASSD_PIF);

figure
X = categorical({'Histogram Matching','Minimum-Maximum','Mean-Standard Deviation','Simple Regression','Pseudo Invariant Features'});
X = reordercats(X,{'Histogram Matching','Minimum-Maximum','Mean-Standard Deviation','Simple Regression','Pseudo Invariant Features'});
Y = [MASSD_HM, MASSD_MM, MASSD_MS, MASSD_SR, MASSD_PIF];
b = bar(X,Y);
xtips = b.XEndPoints;
ytips = double(b.YEndPoints);
labels = string(b.YData);
text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
title('Mean Absolute Spectral Signature Differences')
ylabel('Mean Absolute Reflectance Difference')