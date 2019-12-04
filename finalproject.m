filename_2='NIHMS41007-supplement-4.mov';
videobj_2=VideoReader(filename_2);
numFrames_2=0;
Video_2=cell([],1);
while hasFrame(videobj_2)
F=readFrame(videobj_2);
numFrames_2=numFrames_2+1;
Video_2{numFrames_2}=F;
end
for ii=1:numel(Video_2)
    TotFrames(:,:,:,ii)=Video_2{ii};
end
TotFrames=TotFrames(:,:,1,:);
t_1=reshape(TotFrames,size(TotFrames,1),size(TotFrames,2),size(TotFrames,4));


 for jj=1:size(t_1,3)
img_1=t_1(:,:,jj); 
img_bw=img_1<150;
img_2=bwareafilt(img_bw,[100 2500]);
img_mask(:,:,jj)=img_2;
 end
 
for qq=1:size(t_1,3)-1
stats_t1=regionprops(img_mask(:,:,qq),t_1(:,:,qq),'Centroid','Area','MeanIntensity');
stats_t2=regionprops(img_mask(:,:,qq+1),t_1(:,:,qq+1),'Centroid','Area','MeanIntensity');
xy1=cat(1,stats_t1.Centroid);
a1=cat(1,stats_t1.Area);
mi1=cat(1,stats_t1.MeanIntensity);
tmp=-1*ones(size(a1));
peaks{1}=[xy1,a1,tmp,mi1];

xy2=cat(1,stats_t2.Centroid);
a2=cat(1,stats_t2.Area);
mi2=cat(1,stats_t2.MeanIntensity);
tmp=-1*ones(size(a2));
peaks{2}=[xy2,a2,tmp,mi2];
peaks_matched=MatchFrames(peaks,2,50);

x1=peaks_matched{1}(:,1);
y1=peaks_matched{1}(:,2);
order1=peaks_matched{1}(:,4);
x2=peaks_matched{2}(:,1);
y2=peaks_matched{2}(:,2);
 for pp=1:size(order1)
     if order1(pp)~=-1
    dis(pp)=((x1(pp)-x2(order1(pp)))^2+(y1(pp)-y2(order1(pp)))^2)^(1/2);
     else
     end
 end
 speed(qq)=mean(dis(pp))*(size(t_1,3)/33);
end
qq=1:size(t_1,3)-2;
plot(qq,speed);


BEFORE = speed(1:307);
BULELIGHT = speed(308:645);
AFTER=speed(646:1014);
AFTER_1 = speed(646:804);
AFTER_2 = speed(805:1014);
plot(1:307,BEFORE,'-k')
hold on
plot(308:645,BULELIGHT,'-r')
hold on
plot(646:804,AFTER_1,'-b')
hold on
plot(805:1014,AFTER_2,'-g')

xlabel('Time/(1/30s)');
ylabel('Mean Velocity/s-1');
legend('BEFORE','BULELIGHT','AFTER_1','AFTER_2')
hold off

M_1=mean(BEFORE);
M_2=mean(BULELIGHT);
M_3=mean(AFTER)
M_4=mean(AFTER_2);

[h1,p1]=ttest2(BEFORE,BULELIGHT);
[h2,p2]=ttest2(BULELIGHT,AFTER);
[h3,p3]=ttest2(BEFORE,AFTER);
[h4,p4]=ttest2(BEFORE,AFTER_2);








