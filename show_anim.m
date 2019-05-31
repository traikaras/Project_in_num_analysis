function [] = show_anim( M ,dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

l = length(M(1,:));
figure;

for i = 1:l
    plot(M(:,i),'linewidth',5)
    xlim([1,length(M(:,1))+2])
    ylim([-3.5,3.5])
    pause(dt)
end

