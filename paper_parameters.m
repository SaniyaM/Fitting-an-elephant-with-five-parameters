%creating fourier transform of the image using parameters from the paper
%%
%number of parameters
n=5;

%fourier coefficients
A_x = [0 0 12 0 -14]; 
B_x = [50 18 0 0 0];
A_y = [-60 0 0 0 0];
B_y = [-30 8 -10 0 0];

%variables
x=0; y=0;  
t=1:0.01:100;

%fourier series
for k=1:n  
x = x+ A_x(k)*cos(2*pi*(k)*t) + B_x(k)*sin(2*pi*(k)*t);
y = y+ A_y(k)*cos(2*pi*(k)*t) + B_y(k)*sin(2*pi*(k)*t);
end

x_t = transpose(x);
y_t = transpose(y);
%%
figure(1)
plot(y, -x);