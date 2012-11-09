************************************************************************
*     4-8-stc.f                                                        *
*     Resoluci¢n de un sistema de ecuaciones diferenciales acopladas.  *
*     Sebasti n Torrente Carrillo                                      *
************************************************************************

*En este ejercicio simplemente resolveremos un sistema de ecuaciones
* diferenciales de primer orden (tres ecuaciones concretamente), para
* unos valores iniciales concretos. Al resolverla veremos que el patron
* que obtenemos es caotico, y la m s m¡nima alteraci¢n de las condiciones
* iniciales variaran por completo la soluci¢n. Adem s, en la
* representaci¢n del diagrama de fases de y ì z tendremos un curioso
* patr¢n similar al de una mariposa.

*Para la resoluci¢n de este problema emplearemos el m‚todo de Runge-
* Kutta. Como el sistema de cuaciones ya es de prmer orden no es
* necesario reducirlo. Este problema tambi‚n se podr¡a resolver mediante
* un m‚todo de sucesiones (como el m‚todo de Gauss-Seidel).

* Primero declaro las constantes a usar en el problema.
      parameter (n=3)  !Numero de ecuaciones a resolver.
      dimension y(n),dydx(n),yout(n)
      real r,b,epsilon,x !Variable independiente y constantes del problema
      external derivs,rk4 !Subrutinas empleadas en el problema, la
      !subrutina derivs contiene las funciones a derivar, mientras que
      !la subrutina rk4 contiene el m‚todo de Runge Kutta.
      common/const/r,b,epsilon !Puesto que las csonstantes las empleamos
      !en las subrutinas pasamos las constantes a estas mediante un
      ! common

* Definimos las constantes que intervienen en el problema:
      r=28.0
      b=8.0/3.0
      epsilon=10.0

* Ponemos las condiciones iniciales:
      y(1)=2.0
      y(2)=10.0
      !y(3)=20.0002
      y(3)=5.0
* Damos los valores iniciales y finales de la variable independiente
* x, as¡ como el numero de pasos que daremos en el problema y con esto
* calcularemos cuanto vale el paso h:
      xmin=0.0
      xmax=30.0
      nmax=1000
      h=(xmax-xmin)/real(nmax)

* Inicializo el valor de x:
      x=xmin
* Y abro el archivo donde guardar‚ los datos:
      open(10,file='4-8-stc.dat',status='unknown')

* Inicio el bucle en el que calcular‚ el valor de los puntos de la
* derivada:
        do i=1,nmax+1
* Escribo los valores que tienen en cada momento la variable independiente
* y las funciones. Lo haremos en el documento que vmaos a utilizar y en
* pantalla para ver si los resultados salen con un valor l¢gico:
          write(10,1000) x,y(1),y(2),y(3)
          write(*,1000) x,y(1),y(2),y(3)
* Llamo a derivs para que me calcule los valores de dy(i), y a Runge-
* Kutta para que me calcule los valores de las fucniones:
          call derivs(x,y,dydx)
          call rk4(y,dydx,n,x,h,yout,derivs)
* Redefinimos los valores de x e y(i)
          x=xmin+i*h
          do j=1,n
             y(j)=yout(j)
          end do
      end do
      
* Cerramos el bucle y el documento:
      close(10)
* Escribimos la sentencia que nos dice que el programa ha finalizado sin
* errores:
      write(*,*) '*****FIN DEL PROGRAMA*****'
* Declaramos un formato para presentar los resultados:
 1000 format(4(2x,g15.7))
      stop
      end

************************************************************************
      subroutine derivs(x,y,dy)
************************************************************************
      dimension y(3),dy(3)
      real r,b,epsilon
      common/const/r,b,epsilon
      dy(1)=epsilon*(y(2)-y(1))
      dy(2)=r*y(1)-y(2)-y(3)*y(1)
      dy(3)=y(1)*y(2)-b*y(3)
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
