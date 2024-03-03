clc 
clear all

%basic infomation
L=0.05; %length L=0.05m
b=0.01; %base thickness b=1cm
w=1; %width w=1m
h=15; %heat transfer coefficient h=15W/(m^2*K)
k=180; %thermal conductivity k=180W/(m*K)
T_infinity=25; %T_∞=25℃
T_base=200; %T_base=200℃
sin_theta=0.0995; %angle sin theta
cos_theta=0.9950;


M=101; %node number
det_x=L/(M-1); %nodal spacing del_x

%From energy balance equation to obtain 99 equations 
A=zeros(M-1);
%row1 to row M-2 are based on the energy balance equation
%row1
A(1,1)=-2+2*1*det_x/L-h*det_x^2/(k*L*sin_theta);
A(1,2)=1-(1+1/2)*det_x/L;
%row2 to row M-2 (M-3*M-3)
for m=2:M-2
   A(m,m-1)=1-(m-1/2)*det_x/L;
   A(m,m)=-2+2*m*det_x/L-h*det_x^2/(k*L*sin_theta);
   A(m,m+1)=1-(m+1/2)*det_x/L;
end
%rowM-1 is based on the Boundary Condition
A(M-1,M-2)=k;
A(M-1,M-1)=-k-h*det_x/sin_theta;

%Constant terms
C = zeros(M-1,1); 
%row 1
C(1,1)=-h*det_x^2/(k*L*sin_theta)*T_infinity-T_base*(1-(1-1/2)*det_x/L);
%row2 to rowM-2, the energy balance equation constant term
for m=2:M-2 
   C(m,1)=-h*det_x^2/(k*L*sin_theta)*T_infinity;
end
%row M-1, the boundary condition constant term
C(M-1,1)=-h*det_x/sin_theta*T_infinity;

%node temperatures
T=A \ C;

%the rate of heat transfer from the fin
T_total=0; %node temperature sum except the base and tip nodes
for i=1:M-2
    T_total=T(i)+T_total;
end
Q_fin=h*w*det_x/cos_theta*(T_base+T(M-1)+2*T_total-(2*M-2)*T_infinity);

Q_max=h*2*w*L*(T_base-T_infinity)/cos_theta;
efficiency=Q_fin/Q_max;















