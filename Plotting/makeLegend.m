function han=makeLegend(colors,markerSize,textStr)

figure(); clf;
AR = [240 200];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[200, 200, AR])
set(gcf,'units','points'); hold on;

nLabels = size(colors,1);

cnt = 1;
for ii = nLabels:-1:1
    s = scatter(1,cnt,'o');
    s.MarkerFaceAlpha   = 0.9;
    s.MarkerEdgeAlpha   = 0.9;
    s.SizeData          = markerSize(ii);
    s.MarkerEdgeColor   = colors(ii,:);
    s.MarkerFaceColor   = colors(ii,:);
    
    t=text(3,cnt,textStr{ii},'fontsize',20);
    cnt = cnt+1;
end
xlim([0 10])
ylim([0,nLabels+1])
%axis tight
axis off
han = gca;