************************************************************************
*     1-5-stc                                                          *
*     Obtenci¢n de ceros en una funci¢n mediante los m‚todos de        *
*     bisecci¢n, el m‚todo de Newton-Raphson y el de la secante.       *
*     Sebasti n Torrente Carrillo                                      *
*     A 14-3-2005                                                      *
************************************************************************

C Lo primero que hacemos es comprobar el lugar aproximado en el que se
C encuentran los ceros que queremos calcular:
C El segundo cero se encuentra en torno a x=0.88, y el cuarto cero en
C torno a x=3

C Primero declaremos las variables:

      parameter (nmax=1000) !Numero de iteraciones maximo
C Puesto que uno de los m‚todos es el de Newton Raphson, adem s de la
C funci¢n hemos de declarar su derivada.
      real f,df
      external f,df
C Declaramos los valores de los ceros, para el caso de la bisecci¢n y
C para el m‚todo de la secante necesitaremos dos puntos, de tal modo
C que usaremos la misma pareja de puntos para ambos m‚todos.
      xl1=0.87
      xr1=0.89

      xl2=2.9
      xr2=3.1
      
      n=0
      epsi=1e-4
C Ya solo queda en el programa principal llamar a las subrutinas y
C escribir los resultados.
      call bisect(xl1,xr1,f,nmax,epsi,xm)
      write(*,*) 'Metodo de la bisecci¢n:'
      write(*,*)'Segundo cero encontrado en',xm,' en ',n,' iteraciones.'

      call newtrph(xo,xl1,nmax,epsi,n,xm)
      write(*,*) 'Metodo de Newton-Raphson'
      write(*,*)'Segundo cero encontrado en',xm,' en ',n,' iteraciones.'

      call secante(xo,xl1,nmax,tol,n,xm)
      write(*,*) 'Metodo de la secante'
      write(*,*)'Segundo cero encontrado en',xm,' en ',n,' iteraciones.'
      
      call bisect(xl2,xr2,f,nmax,epsi,xm)
      write(*,*) 'Metodo de la bisecci¢n:'
      write(*,*)'Cuarto cero encontrado en',xm,' en ',n,' iteraciones.'

      call newtrph(xo,xl2,nmax,epsi,n,xm)
      write(*,*) 'Metodo de Newton-Raphson'
      write(*,*)'Cuarto cero encontrado en',xm,' en ',n,' iteraciones.'

      call secante(xl,xl2,nmax,epsi,n,xm)
      write(*,*) 'Metodo de la secante'
      write(*,*)'Cuarto cero encontrado en',xm,' en ',n,' iteraciones.'
      
      write(*,*)
      write(*,*) '********** ESTE PROGRAMA HA FINALIZADO **********'

      stop
      end

************************************************************************
**************************SUBRUTINAS************************************
************************************************************************

***********************Bisecci¢n****************************************
      subroutine bisect(xl,xr,f,nmax,epsi,xm)
      
      fl=f(xl)
      fr=f(xr)
      
*      if (fl*fr .gt. 0.0) then goto 100
      do n=1,nmax
        dif=abs(xl-xr)
        xm=(xr+xl)/2.0
        fm=f(xm)
        if(dif.gt.epsi) then
          if (fm*fr.le.0.0)then
            xl=xm
            fl=fm
          else
            xr=xm
            fr=fm
          endif
        else
          return
        endif
      enddo
      write(*,*)'No se alcanz¢ la precisi¢n deseada.'
      return
*100   write(*,*)'Datos iniciales mal elegidos.'
*      return
      end
************************************************************************

*********************METODO NEWTON-RAPHSON******************************
      subroutine newtrph(xo,x,nmax,epsi,n,x2)

      external f

       do n=1,nmax
         xm=x-f(x)/df(x) !Este es el m‚todo de Newton Raphson
                         !propiamente dicho.
         dif=abs(x-x2)
         if (dif.le.epsi) then
         goto 100
         else
          xo=x
          x=x2
         end if
       end do
100   return
      end
************************************************************************

*********************METODO DE LA SECANTE*******************************
       subroutine secante (xo,x,nmax,epsi,n,x2)

      external f

       do n=1,nmax
         x2=x-((f(x)*(x-xo))/(f(x)-f(xo))) !El m‚todo de la secante
         dif=abs(x-x2)
         if (dif.le.epsi) then
         goto 100
         else
          xo=x
          x=x2
         end if
       end do
100   return
      end
************************************************************************

************************************************************************
****************************FUNCIONES***********************************
************************************************************************
C Declaramos la funcion en la que queremos hallar los ceros.
      function f(x)
      pi=4.0*atan(1.0)
      f=2.0*exp(-2.0*x)-sin(pi*x)
      return
      end
C Necesitaremos la derivada para usar el m‚todo de Newton Raphson
      function df(x)
      pi=4.0*atan(1.0)
      df=-4.0*exp(-2.0*x)-pi*cos(pi*x)
      return
      end

C Como se pude comprobar una vez ejecutado, el valor del cero se
C alcanza con mucha rapidez en todos los m‚todos. Esto se debe
C a la proximidad de los puntos iniciales escogidos.
