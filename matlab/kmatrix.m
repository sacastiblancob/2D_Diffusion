Dx=0.25;
Dy=0.25;
u=0;
v=0;

%Discretizacion espacial
a=4;
b=4;
dx=(8/11);
dy=(8/11);


n=round((a-(-b))/dx+1);
m=round((a-(-b))/dy+1);
x=-a:dx:b;
y=-a:dy:b;
X=kron(ones(1,m),x);
Y1=kron(y,ones(1,n));
Y=Y1(n*m:-1:1);
MX=reshape(X,m,n)';
MY=reshape(Y,m,n)';

%Discretizacion temporal
dt=.025;
tinicial=1+dt;
tfinal=10;
t=tinicial+dt:dt:tfinal;

%Constantes de la ecuacion diferencial computacional
Sx=(Dx*dt)/(dx^2);
Sy=(Dy*dt)/(dy^2);
CFLx=(u*dt)/dx;
CFLy=(v*dt)/dy;
CFL = max(CFLx,CFLy);
%Comprobacion de los parametros Sx, Sy, CFLx y CFLy
%     if CFL <=1
%         disp("Los valores de Sx, Sy, CFLx y CFLy estan dentro de los limites adecuados")
%     else
%         disp("ATENCION: Los valores de Sx, Sy, CFLx y CFLy no estan dentro de los limites adecuados")
%         return
%     end

%Valores de los nodos que estan en la frontera
front=[1:n,n+1:n:n*m-2*n+1,2*n:n:n*m-n,n*m-n+1:n*m];
front=sort(front);
fronts=(n+3:2*n-2);
frontd=(3*n-1:n:n*m-2*n-1);
frontin=(n*m-2*n+3:n*m-n-2);
frontiz=(2*n+2:n:n*m-3*n+2);

%Ensamblaje matriz sistema de ecuaciones para el paso 2 del planteamiento
%numerico
ap=(1+2*Sx+2*Sy);
ax=-Sx;
ay=-Sy;

K11=diag(ones(n-2,1));
K12=ax*diag(ones(n-3,1),1)+ax*diag(ones(n-3,1),-1)+ap*diag(ones(n-2,1));
K21=diag(ones(n-3,1),1)+diag(ones(n-3,1),-1);
K22=ay*diag(ones(n-2,1));
K=kron(K11,K12)+kron(K21,K22);

[Kv,Kr,Kc] = full2csc(K);












