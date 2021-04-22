%% soluting graphing
res = csvread('solution(single).csv');
surf(res);
title('Solution of transport equation');
ylabel('timeSteps (time = 1)'), xlabel('spaceSteps (space = 1)'), zlabel('u(x, t)');