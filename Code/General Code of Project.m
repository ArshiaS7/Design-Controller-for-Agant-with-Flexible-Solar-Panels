%% group #??
%% solving non linear equations for a sample tourque
clear;close all;clc;
% Initializing Parameters
Mq=18.8892;
Mtq=3.6877;
k = 30767.41;
c = 1.525;
It = 4.858;
% Initializing Steps
Time = 100;
dt = 0.00001;
t = 0:dt:Time-dt;
% initializing Tourque
T = zeros(size(t)); %tourque
for n = 1:length(t)/2
    T(n) = 1;
end
for n = length(t)/2+1:length(t)
    T(n) = -1;
end
% initial Conditions for state Variables
x1= zeros(size(t)); % first var of state  = angle
x2= zeros(size(t)); % second var of state = angular velocity
x3= zeros(size(t)); % third....   = q
x4= zeros(size(t)); % fourth...   = q'

% reviewing state space
%  x1' = x2
% x2' = (Mq*T - 2*Mq^2*x2*x3*x4 - Mq*Mtq*x2^2*x3 + k*Mtq*c*x4)/(Mq^2*x3^2 + It*Mq - Mtq^2)
% x3' = x4
% x4' = (Mq*It*x2^2*x3 + Mq^2*x2^2*x3^3 - It*k*x3 - Mq*k*x3^3 - It*c*x4 - Mq*c*x3^2*x4 - Mtq*T + 2*Mq*Mtq*x2*x3*x4)/(Mq^2*x3^2 + Mq*It-Mtq^2)
% simulation
x1_dot = zeros(size(t));
x2_dot = zeros(size(t));
x3_dot = zeros(size(t));
x4_dot = zeros(size(t));
for n = 2:length(t)
    x1_dot(n-1) = x2(n-1);
    x2_dot(n-1) = (Mq * T(n-1) - 2* Mq^2 * x2(n-1) * x3(n-1) * x4(n-1) - Mq * Mtq * x2(n-1)^2 * x3(n-1) + k * Mtq * c * x4(n-1))/(Mq^2 * x3(n-1)^2 + It*Mq - Mtq^2);
    x3_dot(n-1) = x4(n-1);
    x4_dot(n-1) = (Mq * It * x2(n-1)^2 * x3(n-1) + Mq^2 * x2(n-1)^2 * x3(n-1)^3 - It * k * x3(n-1) - Mq * k * x3(n-1)^3 - It * c * x4(n-1) - Mq * c * x3(n-1)^2 * x4(n-1) - Mtq * T(n-1) + 2 * Mq * Mtq * x2(n-1) * x3(n-1) * x4(n-1))/(Mq^2 * x3(n-1)^2 + It*Mq - Mtq^2);
    x1(n) = x1_dot(n-1)*dt + x1(n-1);
    x2(n) = x2_dot(n-1)*dt + x2(n-1);
    x3(n) = x3_dot(n-1)*dt + x3(n-1);
    x4(n) = x4_dot(n-1)*dt + x4(n-1);
end
% plotting
figure(1);
subplot(2,1,1);
plot(t,x2);
title('derivative of angle');
xlabel('time(s)');
ylabel('W(radian/s)');
subplot(2,1,2);
plot(t,x4);
title('derivative of q');
xlabel('time(s)');
ylabel('');

figure(2);
subplot(2,1,1);
plot(t,x1);
title('angle(radian)');
xlabel('time(s)');
ylabel('angle(radian)');
subplot(2,1,2);
plot(t,x3);
title('q');
xlabel('time(s)');
ylabel('q');
%% linearization
clear;clc;close all;
% Initializing Parameters
Mq  = 18.8892;
Mtq = 3.6877;
k   = 30767.41;
c   = 1.525;
It  = 4.858;
% Jacobian for A
syms x1 x2 x3 x4 T;
f1 = x2;
f2 = (Mq*T - 2*Mq^2*x2*x3*x4 - Mq*Mtq*x2^2*x3 + k*Mtq*x3+Mtq*c*x4)/(Mq^2*x3^2 + It*Mq - Mtq^2);
f3 = x4;
f4 = (Mq*It*x2^2*x3 + Mq^2*x2^2*x3^3 - It*k*x3 - Mq*k*x3^3 - It*c*x4 - Mq*c*x3^2*x4 - Mtq*T + 2*Mq*Mtq*x2*x3*x4)/(Mq^2*x3^2 + Mq*It-Mtq^2);

A = jacobian([f1,f2,f3,f4],[x1,x2,x3,x4]);
B = jacobian([f1,f2,f3,f4],T);
A = simplify(A,'Steps',50);
B = simplify(B,'Steps',50);
vpa(A,10)
vpa(B,10)
% Balance Point
A1 = subs(A , [x2,x3,x4],[0,0,0]);
B1 = subs(B,[T , x3],[0 , 0]);
vpa(A1,8)
vpa(B1,8)
%% first root locus and step response for just a simple unity feedback
clear;clc;close all;
s = tf('s');
%entering state space params
A = [0 1 0 0;
    0 0 1451.565 0.0719;
    0 0 0 1;
    0 0 -1912.222 -0.0948]
B = double([0;
    0.242;
    0;
    -0.047])
C = [1 0 0 0];
D = 0;
sys = ss(A,B,C,D);
h = tf(sys)
figure(1)

rlocus(h); %rlocus of h(angle) with just a unity feedback and without any control action
title('\theta')
C = [0 0 1 0];
sys = ss(A,B,C,D);
q = tf(sys)
figure(2);

rlocus(q); %rlocus of q with just a unity feedback and without any control action
title('q')
%% pole placement by state feedback and its root locus
clear;clc;close all;
s = tf('s');
%entering state space params
A = [0 1 0 0;
    0 0 1451.565 0.0719;
    0 0 0 1;
    0 0 -1912.222 -0.0948]
B = double([0;
    0.242;
    0;
    -0.047])
C = [1 0 0 0];
D = 0;
det_controlable= det(ctrb(A,B))
det_observable = det(obsv(A,C))
K = [9.1 9.4 -157.1 -1268.6]; % calculated for desired poles
A_new = A - B*K
det_observable_new = det(obsv(A_new,C))
sys = ss(A_new,B,C,D);
h = tf(sys)
C = [0 0 1 0];
sys = ss(A_new,B,C,D);
q = tf(sys)
rlocus(h); %rlocus after using state feedback
title('\theta')
figure(2);
rlocus(q);
title('q')
% it is obvious that we should check q after using control action to certify its damping action 
%% feedback and control action with state feedback structure(vulnerable to disterbance and Kp!=INF)
clear;clc;close all;
s = tf('s');
%entering state space params
A = [0 1 0 0;
    0 0 1451.565 0.0719;
    0 0 0 1;
    0 0 -1912.222 -0.0948]
B = double([0;
    0.242;
    0;
    -0.047])
C = [1 0 0 0];
D = 0;
K = [9.1 9.4 -157.1 -1268.6]; % calculated for desired poles
A_new = A - B*K
sys = ss(A_new,B,C,D);
h = tf(sys)
C = [0 0 1 0];
sys = ss(A_new,B,C,D);
q = tf(sys)
Gc = (s+2) / (s+25.47);
figure(1);
rlocus(h*Gc); %rlocus after using state feedback

hold on;
ezplot('y=-x',[-50 50]);
title('Root locus of \theta with y=-x'); % K= 730
figure(2);
rlocus(q*Gc); % for K = 730 , its stable.
title('q')
%% protecting from disterbance and enabling system to follow its input
clear;clc;close all;
s = tf('s');
%entering state space params
A = [0 1 0 0;
    0 0 1451.565 0.0719;
    0 0 0 1;
    0 0 -1912.222 -0.0948]
B = double([0;
    0.242;
    0;
    -0.047])
C = [1 0 0 0];
D = 0;
K = [9.1 9.4 -157.1 -1268.6]; % calculated for desired poles
A_new = A - B*K
sys = ss(A_new,B,C,D);
h = tf(sys)
C = [0 0 1 0];
sys = ss(A_new,B,C,D);
q = tf(sys)
Gc = (s+2)^2 / ((s+25.47)*s);
rlocus(h*Gc); %rlocus after using extra control action (s+2)/s % K = 791
hold on;
ezplot('y=-x',[-50 50]);
title('Root locus of H with y=-x'); % K= 791
figure(2);
rlocus(q*Gc); % for K = 791 , its stable.
%% solving non-linear equation with control action and state feedback
%empty