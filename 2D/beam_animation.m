function [] = beam_animation(C,p,U,nt,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fig = figure;
grid on
% x0 = 300;
% y0 = 300;
% width = 1500;
% height = 600;


ptch = patch('faces',C,'vertices',[p(1:2:end,1),p(2:2:end,1)], ...
        'facecolor','interp', 'CData',U(2:2:2*n,1)-p(2:2:2*n-1) ,...
        'edgecolor',[.2,.2,.2]) ;
    


colorbar
%caxis([-0.2,0.6])

%set (gcf, 'position' , [x0, y0, width, height])

for i=2:nt
    if ~ishandle(fig)
        break
    end
    ptch.Vertices = [p(1:2:end,i),p(2:2:end,i)];
    ptch.CData = U(1:2:2*n,i);
%     patch('faces',C(:,1:3),'vertices',reshape(p(:,i),[2,n])', ...
%         'facecolor','w', ...
%         'edgecolor',[.2,.2,.2]) ;
    %hold on; axis image off;
%     patch('faces',edge(:,1:2),'vertices',node, ...
%         'facecolor','w', ...
%         'edgecolor',[.1,.1,.1], ...
%         'linewidth',1.5) ;
    title(['i=' num2str(i)])
    drawnow
    pause(1/2)
    %clf;
    
end

