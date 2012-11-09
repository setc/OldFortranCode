************************************************************************
*     3-8-stc.f                                                        *
*     Calculo de una integral doble hasta la cuarta cifra significativa*
*     Sebasti n Torrente Carrillo                                      *
************************************************************************
C En este ejercicio simplemente emplearemos lo que hemos aprendido a la
C hora de realizar integrales numericas realizando una integral doble.

C Declaramos variables:
      real f,xl,xu,sol!,tol
      common/equis/xx
!      common/tolera/tol
      external f
C Ponemos los l¡mites de integraci¢n de la segunda integral
      xl=0.0
      xu=1.0
!      tol=1e-4
C Realizamos la segunda integral:
      call qgaus(f,xu,xl,sol)
      write(*,*) sol !Escribimos la soluci¢n
      stop
      end
*****************************FUNCIONES**********************************
C La siguiente funci¢n es donde se calcula la primera integraci¢n.
      function f(x)
      common/equis/xx
!      common/tolera/tol
      real fy,yl,yu,sy!,tol
      external fy !Llamamos a la funci¢n a integrar
      yl=x**2  !En esta integraci¢n los l¡mites dependen de x,
      yu=x     !obtendremos como soluci¢n una funci¢n.
!      yl=0.0  !Estas lineas existen como prueba en caso de error.
!      yu=0.1
      call qgaus(fy,yl,yu,sy)
      f=sy  !Esta es la funci¢n que integraremos en el m¢dulo MAIN
      xx=x  !Guardamos el valor de x en la variable xx
      return
      end

C Esta es la funci¢n a integrar, cambiando esta y los l¡mites podemos
C integrar cualquier funci¢n siempre y cuando sea posible.
      function fy(y)
      real x,y
      common/equis/x
      fy=exp(x*y)
      return
      end

*****************************SUBRUTINA**********************************
C Empleamos la subrutina qgaus, gabq es menos precisa, sin embargo hace
C falta declarar m s variables y hay que declarlo todo a doble precisi¢n
C lo que resta legibilidad al programa.
      SUBROUTINE qgaus(func,a,b,ss)
      REAL a,b,ss,func
      EXTERNAL func
      INTEGER j
      REAL dx,xm,xr,w(5),x(5)
      SAVE w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     *.0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666,
     *.9739065285/
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      ss=0
      do 11 j=1,5
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
11    continue
      ss=xr*ss
      return
      END
