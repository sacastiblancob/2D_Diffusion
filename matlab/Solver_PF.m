clear; clc;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  PONTIFICIA UNIVERSIDAD JAVERIANA
%  SERGIO CASTIBLANCO
%  METODOS NUMÉRICOS AVANZADOS
%  SOLVER PARA DIFUSIÓN 2D
%  PARA PROBAR RUTINAS CSC
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Diffusion coefficients
Dx=0.25;
Dy=0.25;
u = 0;
v = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discretizacion espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Limits of the domain
xmin=-4;
xmax=4;
ymin=-4;
ymax=4;

%Number of nodes in X direction (n) and Y direction(m)
n=101;
m =101;

%Space differentials
dx = (xmax-xmin)/(n-1);
dy = (ymax-ymin)/(m-1);

%X and Y vectors
x=xmin:dx:xmax;
y=ymin:dy:ymax;

%X and Y matrices
X=kron(ones(1,m),x);
Y1=kron(y,ones(1,n));
%Y=Y1(n*m:-1:1);
Y = Y1;
MX=reshape(X,n,m)';
MY=reshape(Y,n,m)';

%Time discretization
% dt=.025;
dt=.0128;
tinicial=2;
tfinal=10;
t=tinicial+dt:dt:tfinal;

%Constants for difussion equation
Sx=(Dx*dt)/(dx^2);
Sy=(Dy*dt)/(dy^2);

%Boundary arrays for knowing which nodes are within it
front=[1:n,n+1:n:n*m-2*n+1,2*n:n:n*m-n,n*m-n+1:n*m];
front=sort(front);
frontin=(n+3:2*n-2);
frontd=(3*n-1:n:n*m-2*n-1);
fronts=(n*m-2*n+3:n*m-n-2);
frontiz=(2*n+2:n:n*m-3*n+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing Diffusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%constants
ap=(1+2*Sx+2*Sy);
ax=-Sx;
ay=-Sy;

% % % K full storage (not even in the best use of matlab sparse)
% K11=diag(ones(m-2,1));
% K12=ax*diag(ones(n-3,1),1)+ax*diag(ones(n-3,1),-1)+ap*diag(ones(n-2,1));
% K21=diag(ones(m-3,1),1)+diag(ones(m-3,1),-1);
% K22=ay*diag(ones(n-2,1));
% K=kron(K11,K12)+kron(K21,K22);
% dimen(iv) = length(K);
% %CondK=cond(K);
% %disp('Condicion de la matriz K es: '); disp(CondK);
% K=sparse(K);
% diagk = sparse(diag(diag(K)));

% K CSC storage
[K11v,K11r,K11c] = csc_diag(ones(m-2,1),0);
h1 = [ax*ones(n-2,1) ap*ones(n-2,1) ax*ones(n-2,1)];
d1 = [-1 0 1];
[K12v,K12r,K12c] = csc_diag(h1,d1);
h2 = [ones(m-3,1) ones(m-3,1)];
d2 = [-1 1];
[K21v,K21r,K21c] = csc_diag(h2,d2);
[K22v,K22r,K22c] = csc_diag(ay*ones(n-2,1),0);

[K1v,K1r,K1c] = csc_kron(K11v,K11r,K11c,K12v,K12r,K12c);
[K2v,K2r,K2c] = csc_kron(K21v,K21r,K21c,K22v,K22r,K22c);

[Kv,Kr,Kc] = csc_sum(K1v,K1r,K1c,K2v,K2r,K2c);
% % clearvars K11v K11r K11c K12v K12r K12c K21v K21r K21c K22v K22r K22c K1v K1r K1c K2v K2r K2c
% % dimen(iv) = length(Kc)-1;

% If Jacobi
% %Jacobi
% [Pv,Qv,Qc,Qr] = csc_prejacobi(Kv,Kr,Kc);

% If SOR
% %SOR
% [Pv,Pr,Pc,Qv,Qr,Qc] = csc_preSOR(Kv,Kc,Kr,1.1);

% If LU (really awfuly inefficient)
%
%     %Computing packed LU decomposition for matrix K (Full storage)
%     LU = packlu(K);
%
%     %Computing packed LU decomposition for matrix K (CSC storage)
% %     [LUv,LUr,LUc] = packlu_csc(Kv,Kr,Kc);
%

% SPY matrix K
%spy(K)
%csc_spy(Kv,Kr,Kc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial condition
Co=analitica(X,Y,tinicial,u,-v,Dx,Dy);
conte=1;
zlims = max(Co);

%Error vector
sizet=size(t);
errorPF=zeros(1,sizet(2));

%Caux for handling boundary conditions
Caux=ones(1,m*n-2*m-2*n+4);
trash=zeros(n*m);

%Time loop
times = 0;
for i=tinicial+dt:dt:tfinal
    Ca=analitica(X,Y,i,u,v,Dx,Dy);
    l=1;
    for j=1:m*n
        if find(j==front)
            trash(j)=0;
        elseif find(j==fronts)
            Caux(l) = Co(j) - ay*Ca(j-n);
            l=l+1;
        elseif find(j==frontd)
            Caux(l) = Co(j) - ax*Ca(j+1);
            l=l+1;
        elseif find(j==frontin)
            Caux(l) = Co(j) - ay*Ca(j+n);
            l=l+1;
        elseif find(j==frontiz)
            Caux(l) = Co(j) - ax*Ca(j-1);
            l=l+1;
        elseif j==n+2
            Caux(l) = Co(j) - ay*Ca(j-n)-ax*Ca(j-1);
            l=l+1;
        elseif j==2*n-1
            Caux(l) = Co(j) - ay*Ca(j-n)-ax*Ca(j+1);
            l=l+1;
        elseif j==n*m-2*n+2
            Caux(l) = Co(j) - ay*Ca(j+n)-ax*Ca(j-1);
            l=l+1;
        elseif j==n*m-n-1
            Caux(l) = Co(j) - ay*Ca(j+n)-ax*Ca(j+1);
            l=l+1;
        else
            Caux(l) = Co(j);
            l=l+1;
        end
    end
    
    if i == tinicial+dt
        C = Caux';
    end
    
    tol = 1E-6;
    tic

    % Matlab method
    % C = K\(Caux');

    % Matlab pcg
    % [x,flag,relres,t] = pcg(K,Caux',tol,100,diagk,speye(length(K)),C);

    % Gaussian elimination
    % C = elgauss(K,Caux');

    % Pack-LU
    % C = solpacklu(LU,Caux');

    % Pack-LU CSC
    % C = solpacklu_csc(LUv,LUr,LUc,Caux');

    % %Jacobi CSC
    % [C,t] = csc_jacobi(Kv,Kr,Kc,Pv,Qv,Qc,Qr,Caux',C,100,tol);

    % %SOR CSC
    % [C,t] = csc_SOR(Kv,Kr,Kc,Pv,Pr,Pc,Qv,Qr,Qc,Caux',C,100,tol);

    % % CG CSC
    [C,t] = csc_CG(Kv,Kr,Kc,Caux',C,100,tol);
    %t
    % % DPCG CSC
    % [C,t] = csc_DPCG(Kv,Kr,Kc,Caux',C,100,tol);

    times = times + 1;
    
    CM = reshape(C,n-2,m-2);
    C1 = [pi*ones(1,m-2);CM;pi*ones(1,m-2)];
    C1 = C1(:)';
    C2 = [pi*ones(1,n),C1,pi*ones(1,n)];
    C2(C2==pi)=Ca(C2==pi);

        %MCc=reshape(C,m,n)';
        MCc = reshape(C2,n,m)';
%         MCa=reshape(Co,m,n)'; 
        surf(MX,MY,MCc); axis([-4 4 -4 4 -1 8]); view(2)
        drawnow;
        zlim([0 zlims])
        %pause(dt)
        %errorPF(conte)=norm(C2-Ca);
        %conte=conte+1;
    Co=C2;
end

%toc
%ultimos valores de concentracion calculados reorganizados en matrices

%     MCc=reshape(C2,m,n)';                            %C calculada reordenada en matriz 
%     MCa=reshape(Ca,m,n)';                           %C analitica reordenada en matriz
%     MdifPF=reshape(Ca-C2,m,n)';                        %Diferencia entra C analitica y C calculada reordenada en matriz
% 
%     % plot(t,log10(errorPF))                          %Grafica de la evolucion del error en el tiempo
%     % xlabel('Tiempo (t)')
%     % ylabel('Log10(error)')
% 
%     conv=log10(dt+dx^2);                             %Orden de convergencia
%     faccon=[Sx,Sy,CFLx,CFLy,Sx+Sy,CFLx+CFLy];        %Factores de convergencia

%end