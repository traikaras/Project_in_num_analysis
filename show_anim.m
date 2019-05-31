function [] = show_anim( M )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

l = length(M(1,:));
figure;

for i = 1:l
    plot(M(:,i))
    xlim([-5,length(M(:,1))])
    ylim([-3.5,3.5])
    pause(0.5)
end

