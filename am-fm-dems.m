% AM demodulator model parameters
N = 24;
M = 0.1*N;
F = 0.1*N*1e6;
W = N*1e6;
t_min = 0;
t_max = (2*pi)/F;
t_step = ((2*pi)/F)/100;
t = t_min:t_step:t_max;

% output signals
x = (M+M*cos(F*t)).*cos(W*t);
z = M+M*cos(F*t);

% differentiators to determine time dependence of amplitude
dif1 = gradient(x);
dif2 = gradient(gradient(x));
dif3 = gradient(gradient(gradient(x)));

% energy operators
R0 = (dif1.^2 - x.*dif2);
R1 = (dif2.^2 - dif1.*dif3);

% time dependence of amplitude
ampl = 3.5;
height = 1.2;
temp_h = height;
steps = 15;
y = zeros(steps,101);
for k=1:steps
    y(k,:) = (R0./sqrt(R1).*cos(W*t*1.7)/ampl)+z*temp_h;
    ampl = ampl+1.2;
    temp_h = temp_h-((height-1)/steps);
end

% FM demodulator model parameters
K = 10^8;
A = 0.1*N;
m = 7;
W2 = N*1e6;
O = 0.1*N*1e6;
t_min_2 = 0;
t_max_2 = (2*pi)/O;
t_step_2 = ((2*pi)/O)/100;
t2 = t_min_2:t_step_2:t_max_2;

% output signals
x2 = @(t2)A*cos(W2*t2+m*sin(O*t2));
y2 = (W2+m*O*cos(O*t2))/K;

intx2_2=@(t2)(N/10).*((N/1000+0.001)*sin((N-N/10).*t2)+0.02.*sin((N+N/10).*t2)+0.045.*sin(N.*t2));
intx2_3=@(t2)(N/10).*(0.001*cos((N-N/10).*t2)+0.002.*cos(N.*t2)+0.0008.*cos((N+N/10).*t2));

% differentiators to determine time dependence of amplitude
int1_2=integral(x2,t_min_2,t_max_2);
int2_2=integral(intx2_2,t_min_2,t_max_2);
int3_2=integral(intx2_3,t_min_2,t_max_2);

% energy operators
Q2=(int1_2).^2-x2(t2).*int2_2;
G2=(int2_2).^2-int1_2.*int3_2;
R=sqrt(G2./Q2);

% time dependence of amplitude
z2 = @(t2)((R+m.*O.*cos(O.*t2))./K); 

% AM demodulator signals graph
figure('NumberTitle', 'off', 'Name', 'AM demodulator signals');
plot(x, 'k-');
hold on
plot(z, 'k+');
hold on
for k=1:steps
    plot(y(k,:), 'k-.');
    hold on
end
hold off
legend('x','z','y','Location','north')
xlim([1 101])
set(gca,'FontSize',14)

% FM demodulator signals graph
figure('NumberTitle', 'off', 'Name', 'FM demodulator signals');
subplot(2,1,1)     
plot(t2, x2(t2), 'k-')
legend('x', 'Location', 'south')
xlim([0 2.6*1e-6])
ylim([-2.5 2.5])
set(gca,'FontSize',12)

subplot(2,1,2)      
plot(t2, y2, 'k+')
xlim([0 2.6*1e-6])
hold on
plot(t2, real(z2(t2)), 'k--')
legend('y','z', 'Location','southwest')
set(gca,'FontSize',12)
hold off

% x(t), y(t), z(t) arrays output
fprintf('x(t) AM demodulator:\n\n')
disp(x)

fprintf('\n\ny(t) AM demodulator:\n\n')
disp(y(1,:))

fprintf('\n\nz(t) AM demodulator:\n\n')
disp(z)

fprintf('\n\nx(t) FM demodulator:\n\n')
disp(real(x2(t2)))

fprintf('\n\ny(t) FM demodulator:\n\n')
disp(y2)

fprintf('\n\nz(t) FM demodulator:\n\n')
disp(real(z2(t2)))

; 
% output of the obtained values of maximum absolute and root mean square errors of output signal 
fprintf('\n\nMax. abs. error of AM demodulator output signals:')
a = abs(z - y(1,:));
disp(max(a))

fprintf('RMS error of AM demodulator output signals:')
sub = (z-y(1,:)).^2;
sum_sub = sum(sub);
s = sqrt(sum_sub/100);
disp(s)

fprintf('Max. abs. error of FM demodulator output signals:')
a2 = abs(z2(t2) - y2);
disp(max(a2))

fprintf('RMS error of FM demodulator output signals:')
sub2 = real((z2(t2)-y2).^2);
sum_sub2 = sum(sub2);
s2 = sqrt(sum_sub2/100);
disp(s2)

% calculation of the worst sensitivity of the AM demodulator parameter
fprintf('AM demodulator worst sensitivity:')
y_sens =  ((R0./sqrt(R1))*(0.1/N).*cos(W*t*1.7)/3.5)+z*1.2;
sens = abs(y_sens - y(1,:));
worst_sens = max(a)/max(sens);
disp(worst_sens)