% soluting graphing
res = csvread('solution5.csv');
surf(res);
check3 = csvread('solution3.csv');
check8 = csvread('solution8.csv');
if check3 == check8
    disp("kek\n");
end
title('Solution of transport equation');
ylabel('timeSteps (time = 1)'), xlabel('spaceSteps (space = 1)'), zlabel('u(x, t)');