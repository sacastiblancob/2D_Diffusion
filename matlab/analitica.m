%FUNCION DE LA ECUACIÓN ANALITICA PARA EL CALCULO DE LA CONCENTRACION EN UN
%DOMINIO COMPUTACIONAL 2D PARA UN TIEMPO t
%
%   para los valores de M=1 ; L=1 ; x0=-4 ; y0=4 ; t0=1
%   (para cambiar cualquiera de las anteriores variables, entrar al codigo
%   de la funcion)
%
%   NOTACION C = analitica(x,y,t,u,v,Dx,Dy)
%
%   Donde   x, vector de coordenadas en direccion x
%           y, vector de coordenadas en direccion y
%           t, escalar con el valor del tiempo para el cual se va a
%           calcular la concentracion
%           u, escalar con el valor de la velocidad en x
%           v, escalar con el valor de la velocidad en y
%           Dx, escalar con el valor del coeficiente de difusion en x
%           Dy, escalar con el valor del coeficiente de difusion en y


function C=analitica(x,y,t,u,v,Dx,Dy)
    M=10;
    L=1;
    x0=0;
    y0=0;
    t0=1;
    C=((M/L)/((4*pi*sqrt(Dx*Dy))*(t-t0))).*exp(-(((x-x0-u.*t).^2)*((4*Dx*(t-t0))^-1))-(((y-y0-v.*t).^2)*((4*Dy*(t-t0))^-1)));
end
%end function