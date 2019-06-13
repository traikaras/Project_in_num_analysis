function [] = show_anim( M ,dt)
%SHOWANIM Is a fast way of visualizing time evolution of a vector stored in
%a matrix where each colun is a different time
%   M is the matrix where the column vectors are stored
%   dt is the length of the time step

l = length(M(1,:));
fig = figure;

for i = 1:l
    
    if ~ishandle(fig)
        break
    end
    
    plot(M(:,i),'linewidth',5)
    title(['t = ' num2str(i*dt)])
    xlim([1,length(M(:,1))+2])
    ylim([-1.5,1.5])
    pause(dt)
end

