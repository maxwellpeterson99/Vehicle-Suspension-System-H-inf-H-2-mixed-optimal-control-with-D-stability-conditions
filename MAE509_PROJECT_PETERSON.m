close all; clear all; clc;
opt = sdpsettings('verbose',0,'solver','sedumi');
eps=1*10^-5;

ms=250;
mu=035;
ks=15000;
kt=150000;
bs=450;
bt=1000;

a=.1;
l=1;
Vspre=60;
Vs=Vspre*1000/3600;
t=linspace(0,10,10000);
inp=a/2*(1-cos(2*pi*Vs*t/l));
for i=1:10000
    if t(i) <=1/Vs
        rinp(i)=a/2*(1-cos(2*pi*Vs*t(i)/l));
    else
        rinp(i)=0;
    end
end
maxrinp=max(rinp);

gaindata=zeros(5,4);

%states
A=[0 1 0 0;
   -ks/ms -bs/ms 0 bs/ms;
   0 0 0 1;
   ks/mu bs/mu -kt/mu -(bs+bt)/mu];

%road excitation
B1=maxrinp*[0; 0; -1; bt/mu];

%control input
B2=[0; 1/ms; 0; -1/mu];

%tracked outputs
C1=[1 0 0 0;
    -ks/ms -bs/ms 0 bs/ms];

%feedforwards
D11=[0;0];
D12=[0;0];

sysunc=ss(A,B1,C1,D11);
figure
impulse(sysunc)
gaindata(1,3)=norm(sysunc,inf);
gaindata(1,4)=norm(sysunc,2);

%Hinf optimal
Y1=sdpvar(4,4);
Z1=sdpvar(1,4);
gam1=sdpvar(1);

mat1 = [Y1*A'+A*Y1+Z1'*B2'+B2*Z1 B1 Y1*C1'+Z1'*D12';
      B1' -gam1*eye(1) D11';
      C1*Y1+D12*Z1 D11 -gam1*eye(2)];
  
F1=[];
F1=[F1, Y1>=eps*eye(4)];
F1=[F1, mat1<=-eps*eye(7)];

optimize(F1,gam1,opt);
Y1=value(Y1);
Z1=value(Z1);
K1=Z1*inv(Y1);
gam1=value(gam1);
gaindata(2,1)=gam1;

Acl1=A+B2*K1;
Ccl1=C1+D12*K1;
sys1=ss(Acl1,B1,Ccl1,D11);
gaindata(2,3)=norm(sys1,inf);
gaindata(2,4)=norm(sys1,2);
figure
impulse(sys1)

%H2 optimal
gam2=sdpvar(1);
X2=sdpvar(4,4);
W2=sdpvar(2,2);
Z2=sdpvar(1,4);

mat21 = [A B2]*[X2;Z2] + [X2 Z2']*[A';B2']+ B1*B1';
mat22 = [X2 (C1*X2+D12*Z2)';C1*X2+D12*Z2 W2];
mat23 = trace(W2);

F2 = [];
F2 = [F2, X2>=eps*eye(4)];
F2 = [F2, mat21<=-eps*eye(4)];
F2 = [F2, mat22>=eps*eye(6)];
F2 = [F2, mat23<=gam2];

optimize(F2,gam2,opt);
gam2= value(sqrt(gam2));
gaindata(3,2)=gam2;
K2= value(Z2)*inv(value(X2));

Acl2=A+B2*K2;
Ccl2=C1+D12*K2;
sys2=ss(Acl2,B1,Ccl2,D11);
gaindata(3,3)=norm(sys2,inf);
gaindata(3,4)=norm(sys2,2);
figure
impulse(sys2)

%Hinf-H2 optimal
gam31=sdpvar(1);
gam32=sdpvar(1);
Y3=sdpvar(4,4);
Z3=sdpvar(1,4,'full');
W3=sdpvar(2,2);

mat31  = Y3*A'+A*Y3+B2*Z3+Z3'*B2'+B1*B1';
mat32  = [Y3  (C1*Y3+D12*Z3)';
       C1*Y3+D12*Z3   W3];
mat33  = [Y3*A'+A*Y3+Z3'*B2'+B2*Z3  B1 Y3*C1'+Z3'*D12';
      B1' -eye(1) D11';
      C1*Y3+D12*Z3 D11 -gam31*eye(2)];

F3 = (mat31 <= 0);
F3 = [F3; mat33 <= 0];
F3 = [F3; trace(W3) <= gam32];
F3 = [F3; Y3>= eps*eye(4)];
F3 = [F3; mat32 >= zeros(6)];

obj3=gam31+gam32;
optimize(F3,obj3,opt);
K3= value(Z3)*inv(value(Y3));
gam31=value(sqrt(gam31));
gam32=value(sqrt(gam32));
gaindata(4,1)=gam31;
gaindata(4,2)=gam32;

Acl3=A+B2*K3;
Ccl3=C1+D12*K3;
sys3=ss(Acl3,B1,Ccl3,D11);
gaindata(4,3)=norm(sys3,inf);
gaindata(4,4)=norm(sys3,2);
figure
impulse(sys3)

%Hinf-H2 with D-Stability specifications
gam41=sdpvar(1);
gam42=sdpvar(1);
Y4=sdpvar(4,4);
Z4=sdpvar(1,4,'full');
W4=sdpvar(2,2);

mat41  = Y4*A'+A*Y4+B2*Z4+Z4'*B2'+B1*B1';
mat42  = [Y4  (C1*Y4+D12*Z4)';
       C1*Y4+D12*Z4   W4];
mat43  = [Y4*A'+A*Y4+Z4'*B2'+B2*Z4  B1 Y4*C1'+Z4'*D12';
      B1' -eye(1) D11';
      C1*Y4+D12*Z4 D11 -gam41*eye(2)];
  
ts=0.5;
p_os=0.4;
tr=0.05;

r=(1.8/tr);
alpha=(4.6/ts);
c=(log(p_os)/pi);

%rise time matrix
mat44=[(-r*Y4),(A*Y4+B2*Z4);(A*Y4+B2*Z4)',(-r*Y4)];

%settling time matrix
mat45=A*Y4+B2*Z4+(A*Y4+B2*Z4)'+2*alpha*Y4;

%mpo matrix
mat46=[A*Y4+B2*Z4+(A*Y4+B2*Z4)',c*(A*Y4+B2*Z4-(A*Y4+B2*Z4)');
    c*((A*Y4+B2*Z4)'-(A*Y4+B2*Z4)),A*Y4+B2*Z4+(A*Y4+B2*Z4)'];

F4 = (mat41 <= 0);
F4 = [F4; mat43 <= 0];
F4 = [F4; trace(W4) <= gam42];
F4 = [F4; Y4>= eps*eye(4)];
F4 = [F4; mat42 >= zeros(6)];

%rise time condition
F4=[F4,mat44<=eps*eye(8)];

%settling time condition
F4=[F4,mat45<=eps*eye(4)];

%mpo condition
F4=[F4,mat46<=eps*eye(8)];

obj4=gam41+gam42;
optimize(F4,obj4,opt);
K4= value(Z4)*inv(value(Y4));
gam41=value(sqrt(gam41));
gam42=value(sqrt(gam42));
gaindata(5,1)=gam41;
gaindata(5,2)=gam42;

Acl4=A+B2*K4;
Ccl4=C1+D12*K4;
sys4=ss(Acl4,B1,Ccl4,D11);
gaindata(5,3)=norm(sys4,inf);
gaindata(5,4)=norm(sys4,2);
figure
impulse(sys4)

gaindata