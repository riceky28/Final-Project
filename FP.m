filename1 = 'NIHMS41007-supplement-3.mov';
vidObj = VideoReader(filename1);
numFrames = 0;
Video = cell([],1);
while hasFrame(vidObj)
    F = readFrame(vidObj);
    numFrames = numFrames + 1;
    Video{numFrames} = F ;
end
for ii = 1:numel(Video)
    TotFrames(:,:,:,ii) = Video{ii};
end
TotMask = TotFrames(:,:,1,:);
TotMask = reshape(TotMask,size(TotMask,1),size(TotMask,2),size(TotMask,4));
TotMasks = TotMask;
TotMask(TotMask<=125) = true;
TotMask(TotMask>125) = false;
TotMask = logical(TotMask);
WormsArea = [100 2500];

for jj = 1:size(TotMask,3)
    Masks(:,:,jj) = imfill(TotMask(:,:,jj),'holes');
    Masks(:,:,jj) = bwareafilt(Masks(:,:,jj),WormsArea);
end

for kk = 1:size(TotMask,3)-1
stats_t1 = regionprops(Masks(:,:,kk),TotMasks(:,:,kk),'Centroid','Area','MeanIntensity');
stats_t2 = regionprops(Masks(:,:,kk+1),TotMasks(:,:,kk+1),'Centroid','Area','MeanIntensity');
xy1 = cat(1,stats_t1.Centroid);
a1 = cat(1,stats_t1.Area);
mi1 = cat(1,stats_t1.MeanIntensity);
tmp = -1*ones(size(a1));
peaks{1} = [xy1,a1,tmp,mi1];
xy2 = cat(1,stats_t2.Centroid);
a2 = cat(1,stats_t2.Area);
mi2 = cat(1,stats_t2.MeanIntensity);
tmp = -1*ones(size(a2));
peaks{2} = [xy2,a2,tmp,mi2];
matched = MatchFrames(peaks,2,50);
peaks_matched{kk} = matched;
end

peaks_matched{1,216} = [];
peaks_matched{1,540} = [];
peaks_matched{1,634} = [];
peaks_matched{1,660} = [];
for ll = 1:size(peaks_matched,2)
    if isempty(peaks_matched{1,ll}) ~= 1
        F1 = peaks_matched{1,ll}{1,1};
        F2 = peaks_matched{1,ll}{1,2};
        F1_Centroids = F1(:,1:2);
        F2_Matches = F1(:,4);
        F2_Centroids = F2(:,1:2);
        for mm = 1:size(F2_Matches,2)
            if F2_Matches(mm) ~= -1
                Disp(mm) = (F1_Centroids(1,mm) - F2_Centroids(1,F2_Matches(mm)))^2 + (F1_Centroids(2,mm) - F2_Centroids(2,F2_Matches(mm)))^2;
            else
                Disp(mm) = (F1_Centroids(1,mm) - F1_Centroids(1,mm))^2 + (F1_Centroids(2,mm) - F1_Centroids(2,mm))^2;
            end
        end
        Mean_Velocity(ll) = mean(nonzeros(Disp))/(28/864);
    else Mean_Velocity(ll) = 0;
    end
end
Prior = Mean_Velocity(1:215);
Stimulus = Mean_Velocity(216:536);
After = Mean_Velocity(537:863);
plot(1:215,Prior,'-k')
hold on
plot(216:536,Stimulus,'-r')
plot(537:863,After,'-b')
hold off
set(gca,'YScale','log')
xlabel('Time/(1/30s)')
ylabel('Mean Velocity/(Pixel/Frame)')
ylim([0 3e+06])
legend('Before Stimulus','Stimulus','After Stimulus')

[h1,p1] = ttest2(Prior,Stimulus);
% h1 = 0; p1 = 0.7702;
[h2,p2] = ttest2(Stimulus,After);
% h2 = 0; p2 = 0.6592;
[h3,p3] = ttest2(Prior,After);
% h3 = 0; p3 = 0.9194;

for nn = 1:size(TotMasks,3)
    stats = regionprops(Masks(:,:,nn),TotMasks(:,:,nn),'Centroid');
    xy = round(cat(1,stats.Centroid));
    Centres{nn} = xy;
end

Centroids = false(480,640,864);
se = strel('disk',5);
for oo = 1:size(Centres,2)
    xy_F = Centres{1,oo};
    for pp = 1:size(xy_F)
        Centroids(xy_F(pp,2),xy_F(pp,1),oo) = true;
    end
    Centroids(:,:,oo) = imdilate(Centroids(:,:,oo),se);
end
Centroids = reshape(Centroids,480,640,1,864);
Centroids = uint8(Centroids*255);
cat_Centroids = cat(3,Centroids,zeros(480,640,2,864));
Worms_Centres = TotFrames + cat_Centroids;

video = VideoWriter('WormLocators.avi');
video.FrameRate = 30;
open(video)
for qq = 1:size(Worms_Centres,4)
    Frame = Worms_Centres(:,:,:,qq);
    writeVideo(video,Frame);
end
close(video);