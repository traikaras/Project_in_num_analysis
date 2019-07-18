function [] = show_anim( M ,dt)
%SHOWANIM Is a fast way of visualizing time evolution of a vector stored in
%a matrix where each colun is a different time
%   M is the matrix where the column vectors are stored
%   dt is the length of the time step

l = length(M(1,:));
fig = figure;

for i = l-1:l
    
    if ~ishandle(fig)
        break
    end
    
    plot(M(:,i),'linewidth',5)
    title(['t = ' num2str(i*dt)])
    xlim([1,length(M(:,1))+2])
    ylim([-0.75,0.75])
    pause(dt)
end
% hold on
% plot([0,l],[-1/3,-1/3],'r')
% title('Static Bending Beam')
% xlabel('$x$','Fontsize',25,'Interpreter','latex')
% ylabel('$w$','Fontsize',25,'Interpreter','latex')
% legend('Beam','-1/3','location','best')