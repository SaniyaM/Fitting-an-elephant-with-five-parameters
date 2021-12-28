%implementation of the model from the paper
A = imread('star.png');
% imshow(A)

A_bw = im2bw(A); %image
%A_bw = imresize(A_r, [100 100]);
%invert binary image for tracing boundary
[r, c] = size(A_bw);
for i=1:c
    for j=1:r 
        if A_bw(j,i)==0
            A_bw(j,i)=1;
        else
            A_bw(j,i)=0;
        end
    end
end

dim = size(A_bw);
col = round(dim(2)/2)-90;
row = find(A_bw(:,col), 1);
boundary_0 = bwtraceboundary(A_bw,[210, 163],'S'); %change row, col
xmin = min(boundary_0(:,1));
index = find(boundary_0(:,1)==xmin, 1);
ymin = boundary_0(index,2);

%using xmin and ymin as starting indices
boundary = bwtraceboundary(A_bw, [80, 325], 'S');
L = length(boundary);
 
% figure(1)
% imshow(A_bw);

A_boundary = zeros(L,2);
A_boundary(:,1) = boundary(:,2);
A_boundary(:,2) = boundary(:,1);
% 119
% 158

% 116
% 132

% 90
% 90
freeman_outline = chaincode(A_boundary,false);
% (y,-x)
figure(2)
plot(A_boundary(:,1),-A_boundary(:,2),'g');
% axis([500 850 -350 -75])

K = length(freeman_outline.code);
del_x = zeros(K,1);
del_y = zeros(K,1);
del_t = zeros(K,1);

for i=1:K
    del_x(i) = sign(6-freeman_outline.code(i)).*sign(2-freeman_outline.code(i));
    del_y(i) = sign(4-freeman_outline.code(i)).*sign(freeman_outline.code(i));
    del_t(i) = 1 + ((sqrt(2)-1)/2)*(1 - (-1)^freeman_outline.code(i));
end

%defining x_p = sum up to p del x links
x_p = zeros(K,1);
y_p = zeros(K,1);
t_k = zeros(K,1);

x_p(1) = del_x(1);
y_p(1) = del_y(1);
t_k(1) = del_t(1);

for i=2:K
    x_p(i) = x_p(i-1) + del_x(i);
    y_p(i) = y_p(i-1) + del_y(i);
    t_k(i) = t_k(i-1) + del_t(i);
end

% figure(3)
% % not a perfect outline of the original image  (eg. for a image containing a circle, this doesn't produce a perfect circle)
% plot(y_p,-x_p,'g');   

    
T = t_k(K);
%number of Fourier coefficients
n=100;


%taking derivative of the x and y functions
dx_t_c = 0;
dy_t_c = 0;
dx_t_s = 0;
dy_t_s = 0;

a_n = zeros(n,1);
b_n = zeros(n,1);
c_n = zeros(n,1);
d_n = zeros(n,1);

for k = 1:n
    dx_t_c = (del_x(1)/del_t(1))*(cos(2*pi*k*t_k(1)/T));
    dy_t_c = (del_y(1)/del_t(1))*(cos(2*pi*k*t_k(1)/T));
    dx_t_s = (del_x(1)/del_t(1))*(sin(2*pi*k*t_k(1)/T));
    dy_t_s = (del_y(1)/del_t(1))*(sin(2*pi*k*t_k(1)/T));
    for p=2:K
        dx_t_c = dx_t_c + (del_x(p)/del_t(p))*(cos(2*pi*k*t_k(p)/T) - cos(2*pi*k*t_k(p-1)/T));
        dx_t_s = dx_t_s + (del_x(p)/del_t(p))*(sin(2*pi*k*t_k(p)/T) - sin(2*pi*k*t_k(p-1)/T));
        dy_t_c = dy_t_c + (del_y(p)/del_t(p))*(cos(2*pi*k*t_k(p)/T) - cos(2*pi*k*t_k(p-1)/T));
        dy_t_s = dy_t_s + (del_y(p)/del_t(p))*(sin(2*pi*k*t_k(p)/T) - sin(2*pi*k*t_k(p-1)/T));
    end
    a_n(k) = (T/(2*(k^2)*(pi^2)))*dx_t_c;
    b_n(k) = (T/(2*(k^2)*(pi^2)))*dx_t_s;    
    c_n(k) = (T/(2*(k^2)*(pi^2)))*dy_t_c;
    d_n(k) = (T/(2*(k^2)*(pi^2)))*dy_t_s;
end



%DC Components

eta = zeros(K,1);
eta(1) = 0;
del = zeros(K,1);
del(1) = 0;

for p=2:K
   eta(p) = x_p(p-1) - (del_x(p)/del_t(p))*t_k(p-1);
   del(p) = y_p(p-1) - (del_y(p)/del_t(p))*t_k(p-1);
end

%not normalised
A = (del_x(1)/(2*del_t(1)))*(t_k(1)^2);
C = (del_y(1)/(2*del_t(1)))*(t_k(1)^2);

for p=2:K
    A = A + (del_x(p)/(2*del_t(p)))*(t_k(p)^2 - t_k(p-1)^2) + eta(p)*(t_k(p) - t_k(p-1));
    C = C + (del_y(p)/(2*del_t(p)))*(t_k(p)^2 - t_k(p-1)^2) + del(p)*(t_k(p) - t_k(p-1));
end

A_0 = A/T;
C_0 = C/T;
t=1:0.1:2*K;

% X_f = zeros(K,n);
% Y_f = zeros(K,n);

X_f = 0;
Y_f = 0;

% a_n(2)=-0;
% a_n(3)=-00;
% a_n(4)=25;
% a_n(5)=3;
% c_n(1) = 60;
% c_n(2)=25;
% c_n(3)=15;
% c_n(4)=15;
% c_n(5)=10;
% t_K_2 = t_k + t_k;
%Fourier series
for i=1:n
X_f =  X_f + a_n(i)*cos((2*i*pi*t_k)/T) + b_n(i)*sin((2*i*pi*t_k)/T);
Y_f =  Y_f + c_n(i)*cos((2*i*pi*t_k)/T) + d_n(i)*sin((2*i*pi*t_k)/T);
% figure(4) 
% plot(Y_f, -X_f, 100, 100)
% hold on;
end
% 
X_f = X_f + A_0;
Y_f = Y_f + C_0;

figure(4)
plot(Y_f, -X_f, 100, 100);
