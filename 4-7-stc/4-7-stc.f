************************************************************************
*     4-7-stc                                                          *
*     P‚ndulo de Wilberforce                                           *
*     Sebasti n Torrente Carrillo                                      *
************************************************************************
*En el pendulo de Wilberforce se alternan dos modos de oscilaci¢n, el
* correspondiente a la una oscilaci¢n en el eje de las z, y una
* oscilaci¢n de torsi¢n. Ambas acopladas. Seg£n los resultados te¢ricos
* se tendr¡an que alternar los dos modos de vibraci¢n.
*
* Para la resoluci¢n de estas ecuaciones diferenciales emplearemos el
* m‚todo de Runge Kutta, por lo que antes tendremos que reducir las
* ecuaciones de segundo orden que obtendremos de la Lagrangiana a uno
* de primer orden.

      parameter (n=4)  !Numero de ecuaciones (y variables)
      dimension y(n),dydx(n),yout(n)
      real m,minerc,k,delta,epsilon,x !Variable independiente y ctes.
      external derivs,rk4 !Declaramos las subrutinas que vamos a usar.
      common/const/m,minerc,k,delta,epsilon

* Defino las constantes que intervienen en el problema (en unidades SI)
      m=0.5       !Masa
      minerc=1e-4 !Momento de Inercia
      k=5.0       !Constante El stica
      delta=1e-3  !Constante Torsional
      epsilon=1e-2!Constante de acoplamiento

* Condiciones iniciales para cada una de las funciones
      y(1)=0.01
      y(3)=0.0
      y(2)=0.0
      y(4)=0.0

* Valor inicial y final de la variable independiente, tambi‚n declaramos
* el numero de pasos para poder calcular el valor de cada paso.
      xmin=0.0
      xmax=200.0
      npasos=5000
      h=(xmax-xmin)/real(npasos)

* Inicializamos el valor de x:
      x=xmin
* Abrimos el documento en el que escribiremos los valores que ir 
* tomando la soluci¢n:
      open(10,file='4-7-stc.dat',status='unknown')

* Iniciamos el bucle que nos calcular  los puntos de las funciones:
        do i=1,npasos+1

* Escribimos en el 4-7-stc.dat los valores de la funci¢n
          write(10,*) x,y(1),y(2),y(3),y(4)
          write(*,*) x,y(1),y(2),y(3),y(4)

*Llamamos a las subrutinas:
          call derivs(x,y,dydx)
          call rk4(y,dydx,n,x,h,yout,derivs)

* Damos nuevos valores a x y a y:
          x=xmin+i*h
          do j=1,n
             y(j)=yout(j)
          end do
      end do

      close(10)
      write(*,*) '*****FIN DEL PROGRAMA*****'
      stop
      end

************************************************************************
******************************SUBRUTINAS********************************

******
      subroutine derivs(x,y,dy)
******
      dimension y(4),dy(4)
      real m,minerc,k,delta,epsilon
      common/const/m,minerc,k,delta,epsilon
      dy(1)=y(2)
      dy(2)=-0.5*(2.0*k*y(1)+epsilon*y(3))/m
      dy(3)=y(4)
      dy(4)=-0.5*(2.0*delta*y(3)+epsilon*y(1))/minerc
      return
      end

************************************************************************
      subroutine rk4(y,dydx,n,x,h,yout,derivs)
*     Numerical Recipes (2a ed.): p. 706
************************************************************************
      integer n,nmax
      real h,x,dydx(n),y(n),yout(n)
      external derivs
      parameter (nmax=50)
      integer i
      real h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
      hh=h*0.5
      h6=h/6.
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
      return
      end
