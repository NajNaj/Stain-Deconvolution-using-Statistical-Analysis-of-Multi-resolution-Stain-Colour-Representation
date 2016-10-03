
%% Read images( Training image,Labeling Image, and TestingImage)
%clc, close all;

    %Load Training images 
    Source='./../TrainingSet/1.tif';
    I=imread(Source);
  
    [h, w, ~] = size(I);
    Io = 255;
    NumOfStains=3;

%% show train and test image
 figure,imshow(Source);title('Testing Image');
    
    %% Construct 3*(h*w) image
    B = double(I(:,:,3));
    G = double(I(:,:,2));
    R = double(I(:,:,1));
    RGB = [reshape(R,1,h*w); reshape(G,1,h*w); reshape(B,1,h*w)];

%% Find optical density of images
OD = -log((RGB+1)/Io);    
%%
%L1 norm for the OD matrix
l1 = sum(abs(OD));
%compute transformation matrix elemets
Ele_1=1/sqrt(2);
Ele_2=1/sqrt(6);
Ele_3=sqrt(2/3);

%% Compute transformation matrix
MAT=[Ele_1 -Ele_1 0;-Ele_2 -Ele_2 Ele_3];
%% Perfrorm transformation into the Maxwellian space 
Max_Matrix=zeros(2,h*w);

for i=1:(h*w)
    Max_Matrix(:,i)=MAT*(OD(:,i)./l1(1,i));
end
%% Compute Means for each classified group
[M, W] = fastica(OD,'numOfIC',3);
Ref_Vecs=M;

%     for j=1:NumOfStains
%         Ref_Vecs(j,1)=(Ele_1*M(1,j))-(Ele_2*M(2,j))+(1/3);
%         Ref_Vecs(j,2)=-(Ele_1*M(1,j))-(Ele_2*M(2,j))+(1/3);
%         Ref_Vecs(j,3)=1-Ref_Vecs(j,1)-Ref_Vecs(j,2);
%     end
for z = 1 : 3
   % normalise vector length     
   len=sqrt((Ref_Vecs(1,z)*Ref_Vecs(1,z))+(Ref_Vecs(2,z)*Ref_Vecs(2,z))+ (Ref_Vecs(3,z)*Ref_Vecs(3,z)));
    if (len ~= 0.0)
        Ref_Vecs(1,z)= Ref_Vecs(1,z)/len;
        Ref_Vecs(2,z)= Ref_Vecs(2,z)/len;
        Ref_Vecs(3,z)= Ref_Vecs(3,z)/len;
    end
end
d=pinv(Ref_Vecs)*OD;

%if d(i,j)<0 set d to zero

%%

 
     %% Show results
      H = Io*exp(-Ref_Vecs(:,1)*d(1,:));
      H = reshape(H', h, w, 3);
      H = uint8(H);

      E = Io*exp(-Ref_Vecs(:,2)*d(2,:));
      E = reshape(E', h, w, 3);
      E = uint8(E);

       Bg = Io*exp(-Ref_Vecs(:,3)*d(3,:));
       Bg = reshape(Bg', h, w, 3);
       Bg = uint8(Bg);


%% Show seperated stain for the sample image
figure,
    subplot(2,2,1) 
    imagesc(I);title('Source');
    subplot(2,2,2) 
    imagesc(H);title('Hematoxylin');
    subplot(2,2,3) 
    imagesc(E);title('Eosin');
    subplot(2,2,4) 
    imagesc(Bg);title('Background');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%% Show Density map for each stain
figure,
    subplot(2,2,1)
    imagesc(I);title('Source');
    subplot(2,2,2)
    imagesc(reshape(d(1,:),h,w));title('Density of  H');
    subplot(2,2,3)
    imagesc(reshape(d(2,:),h,w));title('Density of  E');
    subplot(2,2,4)
    imagesc(reshape(d(3,:),h,w));title('Density of  B');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
