% Create data (30k 3D points)
n = 300; 
x = rand(1,n)*100; 
y = rand(1,n)*100; 
z = rand(1,n)*100; 
% Define color of each point based on z value
colors = jet(numel(z));  
[~, ~, depthIdx] = unique(abs(z)); 
colorData = colors(depthIdx,:); 
% Define size gradient based on z value
dotSize = linspace(10,100,numel(x)); 
sizeData = dotSize(depthIdx); 
% Plot 3D points
ax1 = subplot(1,2,1)
scatter3(x,y,z,sizeData,colorData,'filled','o','MarkerEdgeColor','k', ...
    'MarkerFaceAlpha',0.5);
grid on
box on
xlabel('x')
ylabel('y')
zlabel('z')

ax2 = subplot(1,2,2)
scatter3(x,y,z,sizeData,colorData,'filled','o','MarkerEdgeColor','k', ...
    'MarkerFaceAlpha',0.5);
grid on
box on
xlabel('x')
ylabel('y')
zlabel('z')

Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);