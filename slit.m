clear all; close all;clc; format long;
syms y;                                     %Symbolic variable
A=input('Enter slit length: ');
a=A/2;
L=input('Enter the value of L :');
M=input('Enter the value of M :');
n=input('Enter the value of N :');
O=input('Enter the value of O :');
N=input('Enter the number of chebychev polynomials used in the approximation : ');
fprintf('Calculating..... Please wait... ');
h = (2*a/N);    %Step size
x = -a:h:a;
x1 = -a:0.01:a;
T = sym(zeros(N+2,1));  %Chebyshev polynomials of first type
U = sym(zeros(N+1,1));  %Chebyshev polynomials of second type
f1 = zeros(N+1,1);
f2 = f1;
T(1) = 1;                                   
T(2) = y;                                   
U(1) = 1;                                   
U(2) = 2*y;
%Finding Chebyshev polynomials with the help of recursive relations
for j=1:N-1
    T(j+2) = 2*y*T(j+1) - T(j);
    U(j+2) = 2*y*U(j+1) - U(j);
    %U(j+2) = y*U(j+1) + T(j+2);
end
T(N+2) = 2*y*T(N+1) - T(N);
T1 = subs(T,y,y/a); %Finding chebyshev polynomials for -a to a
U1 = subs(U,y,y/a);
I = zeros(N+1,N+1);
A = zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        I(i,j) = int(L*M*n*pi*T1(j+1)*( a/(sqrt(a^2-y^2))),-a,x(i));
        A(i,j) = subs(U1(j),x(i)) - I(i,j);
    end
    f1(i) = L*O*sqrt(a^2 - x(i)^2);     %f(x) = -LO*sqrt(a^2 - x^2)
 
end
c1 = A\f1;    %Gaussian elimination for a system of linear equations
c2 = A\f2;
v1=0;
for i=1:N+1
    v1 = v1 + c1(i)*U1(i);    
end
sol1 = subs(v1,x1); %Finding solution profile by substituting discrete values 
plot(x1,sol1,'r');  %Ploting the slit crack profile 
xlabel('-a < x < +a') 
ylabel('v')
title('Vertical Aperture opening of an uni-axial slit')
xlim([-a-1,a+1]);
ylim([-1,1]);
grid on
hold on
syms x
fprintf('The function V(x) is :');
vpa(simplify(subs(v1,x)),3)