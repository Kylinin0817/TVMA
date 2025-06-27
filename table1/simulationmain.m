clc;
clear;
result1 = zeros(8,9);
result2 = zeros(8,9);
result3 = zeros(8,9);
parpool('local', 16);
for h = 1:4
    [result1(2*h-1:2*h,:)] = simulationfun(100,h,5);
end
for h = 1:4
    [result2(2*h-1:2*h,:)] = simulationfun(300,h,10);
end
for h = 1:4
    [result3(2*h-1:2*h,:)] = simulationfun(500,h,10);
end
delete(gcp('nocreate'));
result1 = result1 ./ result1(:,3);
result2 = result2 ./ result2(:,3);
result3 = result3 ./ result3(:,3);
result1(:,3) = [];
result2(:,3) = [];
result3(:,3) = [];

writematrix(result1, 'outputtable1.xlsx', 'Sheet', 'results1');
writematrix(result2, 'outputtable1.xlsx', 'Sheet', 'results2');
writematrix(result3, 'outputtable1.xlsx', 'Sheet', 'results3');