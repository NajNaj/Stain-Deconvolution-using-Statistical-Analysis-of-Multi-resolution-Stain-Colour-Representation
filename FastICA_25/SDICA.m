
%% add path

%% read source image
Source='./../TrainingSet/6.tif';
Im=imread(Source);

%% show image
figure,
imshow(Im);title('training Image');
axis image

%==========================================================================
% Lifting Dyadic Wavelet Transform
%==========================================================================
[h, g, ht, gt] = SelectLDyWaveletFilterso(2, 4, 1); % h,g for analysis , ht,gt for reconstruction
                                          %1, 2, 1
L =1; % Scale 
[Ls,yh, maxr, maxc] = DecomposeUsingDyadicWTo(double(rgb2gray(Im)),h,g,L);
[ml, nl] = size(Ls);
LL = Ls(maxr +1:ml- maxr, maxc+1:nl- maxc);
HH= Extract_Subbands_HH_rco(L, yh, ml, nl, maxr, maxc);
figure,
 imagesc(HH); title('HH');%% Gray  
figure,
 imagesc(LL); title('LL');%% Gray  
 
