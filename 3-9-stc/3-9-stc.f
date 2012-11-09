************************************************************************
*     3-9a-stc.f                                                       *
* Calculo de los puntos de la funci¢n de distribuci¢n Maxwell-Boltzmann*
*      para diferentes temperaturas                                    *
************************************************************************
C En este programa calculamos los puntos para su represntaci¢n gr fica
C en GNUplot. Para ello definiremos la funci¢n dependiente de dos
C variables (T,V) y calcularemos los puntos para canda temperatura.

C Primero declaramos y definimos las variables que necesitamos.
      real f,T1,T2,T3,h,v,v0,vmax
      external f !Declaramos la funci¢n
      integer i,imax
      v0=-100.0   !Hemos hecho una prueba previa con los valores en los
      vmax=100.0  !que evaluamos la funci¢n para encontrarn el m s
      imax=1000  !apropiado.
      T1=100.0
      T2=300.0
      T3=500.0
      h=(vmax-v0)/imax !Calculamos el paso.

C Ahora generamos tres archivos para guardar los datos necesarios para
C representar graficamente cada funci¢n, a cada archivo *.dat se le da
C el nombre del ejercicio con la letra a,b o c, correspondiendo las
C temperaturas T1=100, T2=300, T3=500 respectivamente.

      open(10,file='3-9a-stc.dat',status='unknown')
C Iniciamos el bucle
      do i=1,imax
         v=v0+h*i !Paso de la variable.
         write(10,*) v,f(T1,v)
         write(*,*) v,f(T1,v)  !Esta linea se usa para comprobar si
                                !el programa va bien.
      enddo
      close(10)
      
      open(11,file='3-9b-stc.dat',status='unknown')
      do i=1,imax
         v=v0+h*i
         write(11,*) v,f(T2,v)
      enddo
      close(11)
      
      open(12,file='3-9c-stc.dat',status='unknown')
      do i=1,imax
         v=v0+h*i
         write(12,*) v,f(T3,v)
      enddo
      close(12)

!Imprimimos en pantalla un mensaje que nos indique el fin de este.
      write(*,*)'***************FIN DEL PROGRAMA****************'
      
      stop
      end
*****************************FUNCIONES**********************************
C A pesar de que la funci¢n depende de dos variables en nuestro caso
C una de las dos (T) tend  unos valores fijos determinados, siendo
C la variable v la £nica variable independiente.
      function f(T,v)
      real pi,m,k  !Definimos las constantes
!      pi=4.0*atan(1.0)
!      m=39.984
      k=1.380657e-23
      f=(1.0/T)**(1.5)*v**2*exp(-(v**2)/T) !Escribimos la funci¢n
      return
      end
