************************************************************************
*     0-12-stc                                                         *
*     Representar conjuntamente las curvas de la ley de Wien y la ley  *
*     de Plank para la intensidad de la radiaci�n de un cuerpo negro.  *
*     Sebasti�n Torrente Carrillo                                      *
*     A 2/3/2005                                                       *
************************************************************************

!Primero definimos todas las constantes que necesitaremos:
      implicit real*8 (a-h,o-z)
      pi=4*atan(1.0)
      c=299792458.0
      h=6.62606876e-34 !He decidido mantener en nombre k a pesar de ser
                       !un caracter reservado para entero.
      k=1.380650e-23

      T1=2.7d0
      T2=300d0
      T3=58000d0
      x=0.0
*      integer i,j
      imax=1000
*Creamos ahora los documentos en los que el programa escribir� los
*datos que ha obtenido.Repetiremos este proceso para cada temperatura.
      open (10, file='0-12-stc1.plt',status='unknown')
        write(10,*) '#lamda Ew(lamd,T=',T1,')Ep(lamd,T=',T1,')'
          do i=1,imax
             real_i=real(i)*8e-6
             write(10,400) real_i,wien(real_i,T1),plck(real_i,T1)
             write(*,*) wien(real_i,T1)
          enddo
        write(10,*)'#Este documento representa los valores de las'
        write(10,*)'#funciones de la ley de Wien y de Planck'
        write(10,*)'#de la radiaci�n de un cuerpo negro en funci�n'
        write(10,*)'#de la longitud de onda lambda para una'
        write(10,*)'#temperatura T=',T1,'K.'
      close (10)

      open (20, file='0-12-stc2.plt',status='unknown')
        write(20,*) '#lamda Ew(lamd,T=',T1,')Ep(lamd,T=',T1,')'
        do i=1,imax
             real_i=real(i)*4e-8
             write(20,400) real_i,wien(real_i,T2),plck(real_i,T2)
          enddo
        write(20,*)'#Este documento representa los valores de las'
        write(20,*)'#funciones de la ley de Wien y de Planck'
        write(20,*)'#de la radiaci�n de un cuerpo negro en funci�n'
        write(20,*)'#de la longitud de onda lambda para una'
        write(20,*)'#temperatura T=',T2,'K.'
      close (20)

      open (30, file='0-12-stc3.plt',status='unknown')
        write(30,*) '#lamda Ew(lamd,T=',T1,')Ep(lamd,T=',T1,')'
          do i=1,imax !As� obtenemos los cien pasos de 0 a diez.
             real_i=real(i)*2e-9
             write(30,400) real_i,wien(real_i,T3),plck(real_i,T3)
          enddo
        write(30,*)'#Este documento representa los valores de las'
        write(30,*)'#funciones de la ley de Wien y de Planck'
        write(30,*)'#de la radiaci�n de un cuerpo negro en funci�n'
        write(30,*)'#de la longitud de onda lambda para una'
        write(30,*)'#temperatura T=',T3,'K.'
      close (30)

      write(*,*)'Este programa ha finalizado.'
      

400   format(f8.5,6x,f10.8,10x,f10.8)
500   format(f10.9,6x,f10.1,10x,f10.1)
600   format(f10.9,6x,f10.1,10x,f10.1)

      stop
      end

      function wien(lamda,temp)
!      wien=(2.0*pi*h*c**2.0)/(lamda**5.0*exp(c*h/(lamda*k*temp)))
       wien=(2.0*pi*h*c**2.0)/lamda**5.0
      return
      end

      function plck(lamda,temp)
!      plck=(2.0*pi*h*c**2.0)/(lamda**5.0*(exp(c*h/(lamda*k*temp)-1.0)))
       plck=1.0/(exp(c*h/(lamda*k*temp)))
      return
      end
