%SOLVER PARA UN DOMINIO COMPUTACIONAL 2D DE LA ECUACION DIFERENCIAL DE
%TRANSPORTE - Metodo pasos fraccionados (PF)
%
%DOMINIO ESPACIAL   -a < x < b
%                   -a < y < b
%
%DOMINIO TEMPORAL   tincial < t < tfinal
%   (donde t inicial es ligeramente mayor a t0)
%
%   En este caso particular se toma t0 = 1s (ver funcion analitica), y
%   tinicial = 1s + dt
%

Dx=0.25;
Dy=0.25;
u=0.1;
v=0.1;

%Discretizacion espacial
a=4;
b=4;
dx=1;
dy=1;
n=(a-(-b))/dx+1;
m=(a-(-b))/dy+1;
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

%Comprobacion de los parametros Sx, Sy, CFLx y CFLy
if (Sx + Sy) <=0.5 && ((CFLx^2/Sx) + (CFLy^2/Sy)) <=2
    disp("Los valores de Sx, Sy, CFLx y CFLy estan dentro de los limites adecuados")
else
    disp("ATENCION: Los valores de Sx, Sy, CFLx y CFLy no estan dentro de los limites adecuados")
    return
end

%Valores de los nodos que estan en la frontera
front=zeros(1,2*n+2*m-4);
cont=n+1;
for i=1:m
    if i==1
        front(1:n)=1:n;
    elseif i==m
        front(2*n+2*m-4-n:2*n+2*m-4)=n*m-n:m*n;
    else
        front(cont)=(i-1)*n+1;
        front(cont+1)=i*n;
        cont=cont+2;
    end
end

%Ensamblaje matriz sistema de ecuaciones para el paso 2 del planteamiento
%numerico

A=zeros(n*m,n*m);
for i=1:n*m
    if find(i==front)
        A(i,i)=1;
    else
        A(i,i)=(1+2*Sx+2*Sy);
        A(i,i-1)=-Sx;
        A(i,i+1)=-Sx;
        A(i,i-n)=-Sy;
        A(i,i+n)=-Sy;
    end
end
CondA=cond(A);
disp('Condicion de la matriz A es: '); disp(CondA);
% A = sparse(A);
% 
% % spy(A)
% 
% %Solver
% 
% Co=analitica(X,Y,tinicial,u,-v,Dx,Dy);
% conte=1;
% sizet=size(t);
% errorPF=zeros(1,sizet(2));
% Caux=zeros(1,m*n);
% tic
% for i=tinicial+dt:dt:tfinal
%     Ca=analitica(X,Y,i,u,-v,Dx,Dy);
%     for j=1:m*n
%         if find(j==front)
%             Caux(j)=analitica(X(j),Y(j),i,u,-v,Dx,Dy);
%         else
%             Caux(j)=Co(j-1)*CFLx + Co(j-n)*CFLy + Co(j)*(1-CFLx-CFLy);      %Paso 1, calculo de Caux explicito
%         end
%     end
%     C = A\(Caux');                                                          %Paso 2, calculo de C(t+1) implicito
%     C=C';
%     
%     %MCc=reshape(C,m,n)';
%     %MCa=reshape(Ca,m,n)'; 
%     %surf(MX,MY,MCa-MCc); axis([-4 4 -4 4 -1 8]); %view(2)
%     %drawnow;
%     %pause(dt)
%     errorPF(conte)=norm(C-Ca);
%     conte=conte+1;
%     Co=C;
% end
% toc
% %ultimos valores de concentracion calculados reorganizados en matrices
% 
% MCc=reshape(C,m,n)';                            %C calculada reordenada en matriz 
% MCa=reshape(Ca,m,n)';                           %C analitica reordenada en matriz
% MdifPF=reshape(Ca-C,m,n)';                        %Diferencia entra C analitica y C calculada reordenada en matriz
% 
% plot(t,log10(errorPF))                          %Grafica de la evolucion del error en el tiempo
% xlabel('Tiempo (t)')
% ylabel('Log10(error)')
% 
% conv=log10(dt+dx^2);                             %Orden de convergencia
% faccon=[Sx,Sy,CFLx,CFLy,Sx+Sy,CFLx+CFLy];        %Factores de convergencia
% 
% %end
