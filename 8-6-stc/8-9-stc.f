************************************************************************
      subroutine ajlineal(a,err_a,b,err_b,ndat)
************************************************************************
!Esta Subrutina genera un ajuste lineal.
      implicit real*8 (a-h,o-z)
*Declaramos los valores de x e y como arrays, como no sabemos a priori
*la cantidad de pares de datos de los que disponemos haremos que sean
*arrays de dimensi¢n n y luego definiremos n al leerlo en el documento.
      dimension x(ndat),y(ndat)
      integer i


*Necesitamos ahora abrir el documento donde tenemos los datos, este
*tiene el nombre '0-11-stc.in' y ha de tener el siguiente formato:
*     n =numero de datos
*     x(1) y(1)
*     x(2) y(2)
*     ...  ...
*     x(n) y(n)
      open(100, 'aajustar.dat', blank='null') !Se¤alamos que los
                                              !espacios no cuentan.
      read(100,*) ndat !Leemos el numero de datos.
*Ahora leeremos los datos de cada array, para ellos empleamos un bucle.

      do i=1,ndat
        read(100,*) x(i),y(i)
      enddo
*Ahora que tenemos los datos podemos obtener a y b con sus errores.
*Como no necesitamos m s el archivo con los datos de entrada lo cerramos
      close 100

*Para calcular los datos necesitamos las medias de x e y.
      a=0.0
      b=0.0
      err_a=0.0
      err_b=0.0
      xmed=0.0
      ymed=0.0
      D=0 !Esta variable la hemos declarado con mayusculas para que tenga
          !el mismo nombre que la D del enunciado.

      do i=1,ndat
        xmed=xmed+x(i)
        ymed=ymed+y(i)
      enddo

      xmed=xmed/real(ndat) !La divisi¢n la realizamos fuera del bucle
      ymed=ymed/real(ndat) !para ahorrarle calculos al ordenador.

*Para calcular a y b ya solo necesitamos D=sum(x(i)-<x>)**2

      do i=1,ndat
        D=D+(x(i)-xmed)**2.0
      enddo

!Ahora calculemos a:

      do i=1,ndat
        a=a+(x(i)-xmed)*y(i) !Primero sumamos todos los t‚rminos
      enddo
      a=a/D                  !Y despu‚s lo dividimos por D.

!Una vez que hemos obtenido a, obtener b es muy facil.

      b=ymed-a*xmed

!Ya solo nos queda calcular los errores err_a y err_b, como ambos
!errores tienen un termino com£n sum[d(i)]/(n-2). POr lo que antes de
!Calcular los errores de a y b calculamos ese sumatorio.

      real d=0
      do i=1,ndat
        d=d+(y(1)-a(x(i))-b)**2.0
      enddo
      d=d/(real(ndat)-2.0)

      err_a=sqrt(d/D)

      err_b=sqrt(d*(xmed**2.0/D+1.0/real(ndat))

      return
      end
