if(idx_1==1)Raster=RasterFS0; s='I cell index';end
if(idx_1==0)Raster=Raster_P0; s='E cell index';end
%%subplot(2,1,1);
for i=1:size(Raster,1)
for spike=1:size(Raster{i},1)
x=Raster{i}(spike);
%%plot([x,x+0.3]/1000,[i,i+1.1],'Color','k','LineWidth',0.2); %%here we plot the result against
plot([x,x]/1000,[i,i],'.','MarkerSize',4,'Color','k'); %%here we plot the result against
%%cell member by using little lines; time (sec)
end
set(gca,'FontSize',13);
ylabel(s,'FontSize',14);
xlim([10 16]); %horizontal range
ylim([-2 size(Raster,1)+5])% here we set the highest limit for the rasterplot
hold on;
end