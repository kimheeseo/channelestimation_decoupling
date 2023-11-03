% No_reflector
N = 10;
x = zeros(N, 4000);

for i=1:N
     fid =fopen(['.\No_reflector\',num2str(i),'.txt'],'r');
     temp = fscanf(fid, '%f\n');
     x(i,:) = temp.^2;
     fclose(fid);
end

% average
y = sum(x)/N;
t = [0:1/4000:1-1/4000];
figure;
plot(t,10*log10(y));
grid on;
% ylim([0 40]);
xlabel('time (us)');
ylabel('Relative power (dB)');
title('Average power delay profile with no reflector');

% with_reflector
N = 20;
x = zeros(N, 2000);

for i=1:N
     fid =fopen(['.\with_reflector\',num2str(i),'.txt'],'r');
     temp = fscanf(fid, '%d\n');
     x(i,:) = (temp-mean(temp(1000:2000))).^2;
     fclose(fid);
end

% average
y = sum(x)/N;
t = [0:1/2000:1-1/2000];
figure;
plot(t,10*log10(y));
grid on;
ylim([0 40]);
xlabel('time (us)');
ylabel('Relative power (dB)');
title('Average power delay profile with one reflector');
