function han = PolarPlot(th,rho,colors,markerSize)
% make polar scatter plot for multiple categories

nCategories     = size(th,2);

z       = rho.*exp(1i*th);
reZ     = real(z);
imZ     = imag(z);
zM      = nanmean(z);

if  nargin<3
    colors = unifrnd(0.5,1,[nCategories 3]);
end
if nargin<4
    markerSize =ones(1,nCategories)*100;
end


han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])

xLims = [-1 1]*max(rho(:));
xlim(xLims);
ylim(xLims); hold on;

lCol = [0.9 0.9 0.9];
plot(xlim,[0 0],'color',lCol,'linewidth',2);
plot([0 0],ylim,'color',lCol,'linewidth',2);
plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);

% reference circle for magnitude
t = 0:0.01:2*pi;
xx = xLims(1)*cos(t)*0.75;
yy = xLims(1)*sin(t)*0.75;
plot(xx,yy,'-','color',lCol,'linewidth',1);
xpos = xLims(2)*0.75; xstr = round(xpos*100)/100;
t=text(xpos,0,num2str(xstr));
t.FontSize=18;
t.HorizontalAlignment='center';
t.VerticalAlignment='top';

for jj = 1:nCategories
    s=scatter(reZ(:,jj),imZ(:,jj),'o');
    s.MarkerFaceAlpha   = 0.7;
    s.MarkerEdgeAlpha   = 0.7;
    s.SizeData          = markerSize(:,jj);
    s.MarkerEdgeColor   = colors(jj,:);
    s.MarkerFaceColor   = colors(jj,:);

end
for jj = 1:nCategories
    plot([0 real(zM(jj))],[0 imag(zM(jj))],'linewidth',5,'color',colors(jj,:)*0.8)
end
axis off

return 