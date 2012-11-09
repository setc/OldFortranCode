************************************************************************
*     3-9b-stc                                                         *
*     Calculo de la velocidad cuadr tica media en la distribuci¢n de   *
*     Maxwell-Boltzman                                                 *
*     Sebasti n Torrente Carrillo                                      *
************************************************************************

C Calcularemos la velocidad cuadr tica media para las tres temperaturas
C para calcularla tenemos antes que crear la funci¢n que corresponde a
C la velocidad cuadr tica media.

C Primero declaramos las variables y las funciones que vamos a emplear.
      real xl,xu,sol,vcT1,vcT2,vcT3
      external vcT1,vcT2,vcT3
      dimension sol(3)
C La integral tiene un intervalo de infinito a infinito, pero al ser
C convergente (si no el ejercicio no tendr¡a sentido) hay un punto en
C que la contribuci¢n de la funci¢n es despreciable, por lo que si
C tomamos unos l¡mites de integraci¢n lo suficientemente amplios como
C para que lleguemos a terminos despreciables.
      xl=-100.0
      xu=100.0

C Llamamos a la subrutina para cada temperatura e imprimimos los
C resultados en pantalla
      call qgaus(vcT1,xl,xu,sol(1))
      write(*,*) '<v**2(T1)>=',sol(1)

      call qgaus(vcT2,xl,xu,sol(2))
      write(*,*) '<v**2(T2)>='sol(2)
      
      call qgaus(vcT3,xl,xu,sol(3))
      write(*,*) '<v**2(T3)>='sol(3)

      stop
      end
      
********************************FUNCIONES*******************************

C Primero definimos la funci¢n correspondiente a la funci¢n de
C distribuci¢n de velocidades
      function f(T,v)
      real pi!,m,k  !Definimos las constantes
      pi=4.0*atan(1.0)
!      m=1.0
!      k=1.0
!      f=pi*(1.0/T)**(3.0/2.0)*v**2*exp(-(v**2/T))
       f=(1.0/T)**1.5*v**2*exp(-(v**2/T))
      return
      end
C Y ahora escribimos las expresiones a integrar. Estas funciones
C dependen de la anterior y de la variable v. La declaramos para cada
C velocidad.
      function vcT1(v)
      external f
      vcT1=v**2*f(100.0,v)**2
      return
      end

      function vcT2(v)
      external f
      vcT2=v**2*f(300.0,v)**2
      return
      end

      function vcT3(v)
      external f
      vcT3=v**2*f(500.0,v)**2
      return
      end

******************************SUBRUTINAS********************************
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
