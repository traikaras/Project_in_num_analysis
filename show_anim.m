function [] = show_anim( M ,dt)
%SHOWANIM Is a fast way of visualizing time evolution of a vector stored in
%a matrix where each colun is a different time
%   M is the matrix where the column vectors are stored
%   dt is the length of the time step

l = length(M(1,:));
figure;

for i = 1:l
    plot(M(:,i),'linewidth',5)
    title(['t= ' num2str(i*dt)])
    xlim([1,length(M(:,1))+2])
    ylim([-3.5,3.5])
    pause(dt)
end

