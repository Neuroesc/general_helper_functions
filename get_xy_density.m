function d = get_xy_density(x,y,r1,r2,k)





[f,xe,ye] = histcounts2(x(:),y(:),[r1 r2]);
fs = imgaussfilt(f,k,'Padding',0,'FilterDomain','spatial');

xg = movmean(xe,2,'Endpoints','discard');
yg = movmean(ye,2,'Endpoints','discard');

[X,Y] = meshgrid(yg,xg);
d = interp2(X,Y,fs,y(:),x(:),'nearest');


% figure
% imagesc(fs)
% 
% 
% keyboard





































