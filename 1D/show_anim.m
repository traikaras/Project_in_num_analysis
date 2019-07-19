function [] = show_anim( M ,dt,mov)
%SHOWANIM Is a fast way of visualizing time evolution of a vector stored in
%a matrix where each colun is a different time
%   M is the matrix where the column vectors are stored
%   dt is the length of the time step

l = length(M(1,:));
fig = figure;
if mov
    currentFolder=pwd;
    writerObj = VideoWriter(strcat(currentFolder,'/beam_stationary.avi'));
    set(writerObj,'FrameRate',1/dt) % Frame rate fixed at 4 as any faster and you miss the first frame
    open(writerObj);
end

for i = 1:l
    
    if ~ishandle(fig)
        break
    end
    
    plot(M(:,i),'linewidth',5)
    title(['t = ' num2str(i*dt)])
    xlim([1,length(M(:,1))+2])
    ylim([-0.75,0.75])
    if mov
        writeVideo(writerObj,getframe(fig));
        if i ==1
            axis tight manual
            set(gca,'NextPlot','replaceChildren')
        end
    
    else
        pause(dt)
    end
  

end
if mov
    close(writerObj);
end
end