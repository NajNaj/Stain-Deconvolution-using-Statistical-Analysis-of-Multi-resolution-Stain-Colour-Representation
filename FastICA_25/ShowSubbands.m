function ShowSubbands(LL,LH,HL,HH,NumLevel)

%min =val , max=val;
         
         LL(:,:,1) = (LL(:,:,1)-min(reshape(LL(:,:,1),1,[])))/max(reshape(LL(:,:,1),1,[])); 
         LL(:,:,2) = (LL(:,:,2)-min(reshape(LL(:,:,2),1,[])))/max(reshape(LL(:,:,2),1,[])); 
         LL(:,:,3) = (LL(:,:,3)-min(reshape(LL(:,:,3),1,[])))/max(reshape(LL(:,:,3),1,[]));

         LH(:,:,1) = (LH(:,:,1)-min(reshape(LH(:,:,1),1,[])))/max(reshape(LH(:,:,1),1,[])); 
         LH(:,:,2) = (LH(:,:,2)-min(reshape(LH(:,:,2),1,[])))/max(reshape(LH(:,:,2),1,[])); 
         LH(:,:,3) = (LH(:,:,3)-min(reshape(LH(:,:,3),1,[])))/max(reshape(LH(:,:,3),1,[]));

         HL(:,:,1) = (HL(:,:,1)-min(reshape(HL(:,:,1),1,[])))/max(reshape(HL(:,:,1),1,[])); 
         HL(:,:,2) = (HL(:,:,2)-min(reshape(HL(:,:,2),1,[])))/max(reshape(HL(:,:,2),1,[])); 
         HL(:,:,3) = (HL(:,:,3)-min(reshape(HL(:,:,3),1,[])))/max(reshape(HL(:,:,3),1,[]));

         HH(:,:,1) = (HH(:,:,1)-min(reshape(HH(:,:,1),1,[])))/max(reshape(HH(:,:,1),1,[])); 
         HH(:,:,2) = (HH(:,:,2)-min(reshape(HH(:,:,2),1,[])))/max(reshape(HH(:,:,2),1,[])); 
         HH(:,:,3) = (HH(:,:,3)-min(reshape(HH(:,:,3),1,[])))/max(reshape(HH(:,:,3),1,[]));
% LL=(LL-min(LL(:)))/max(LL(:));
% LH=(LH-min(HL(:)))/max(HL(:));
% HL=(HL-min(HL(:)))/max(HL(:));
% HH=(HH-min(HH(:)))/max(HH(:));



%% show Bands
figure, 
h(1)=subplot(2,2,1);imagesc(LL);str= sprintf('Approximation Level:%d',NumLevel);title(str);
colorbar;axis image
h(2)=subplot(2,2,2);imagesc(LH);str= sprintf('Horizantal Level:%d',NumLevel);title(str);
colorbar;axis image
h(3)=subplot(2,2,3);imagesc(HL),str= sprintf('Vertical Level:%d',NumLevel);title(str);
colorbar;axis image
h(4)=subplot(2,2,4);imagesc(HH),str= sprintf('Diagonal Level:%d',NumLevel);title(str);
colorbar;axis image
linkaxes(h,'xy');
set(gcf,'units','normalized','Position', [0 0 1 1]),
%print(sprintf('/Users/najahsubaie/Desktop/PhD_2015/MaxellResults/%itif',NumLevel),'-dpng','-r0');

% F = getframe(gcf);
% imwrite(F.cdata, './SBICA_Images/A1.png');
% close all;














end