************************************************************************
*     2-6-stc                                                          *
*     Calculo de una resistencia equivalente.                          *
*     Sebasti n Torrente Carrillo                                      *
*     31-3-2005                                                        *
************************************************************************

C En este programa calcularemos la resistencia equivalente de una red
C de resistencias dadas empleando las reglas de Kirchoff y los metodos
C de resoluci¢n de ecuaciones lineales.

C El primer paso corresponde a el calculo de las ecuaciones que
C obtengamos empleando las reglas de Kirchoff y su representaci¢n
C en forma matricial.

C Declaramos la dimensi¢n de la matriz juntoa otras variables.
      integer i,j,np,n  ! i y j son indices para los bucles.
      parameter (n=11,np=11)
      dimension r(n-1),y(n,n),b(n),b0(n),indx(np)
      real v

C Damos valores a v y a los distintos tipos de resistencia, los
C valores que pueden tomar las resistencias se explican en el
C documento adjunto:

      v=100.0 !Valor arbitrario.
      res1=50.0
      res2=70.71078
      res3=100.0

C Inicializamos todos los valores de la matriz a cero
      do i=1,n
        b(i)=0.0
        r(i)=0.0
        do j=1,n
           y(i,j)=0.0
        enddo
      enddo

C Y ahora asignamos valores:
      b(11)=v

      r(1)=res1
      r(2)=res1
      r(3)=res3
      r(4)=res2
      r(5)=res2
      r(6)=res2
      r(7)=res2
      r(8)=res2
      r(9)=res1
      r(10)=res1
C Definimos los valores de la matriz distintos de cero:
      y(1,1)=1.0
      y(1,2)=1.0
      y(1,11)=-1.0

      y(2,1)=-1.0
      y(2,3)=1.0
      y(2,4)=1.0
      y(2,6)=1.0

      y(3,2)=-1.0
      y(3,4)=-1.0
      y(3,5)=-1.0
      y(3,7)=-1.0

      y(4,5)=-1.0
      y(4,7)=1.0
      y(4,8)=-1.0
      y(4,10)=-1.0

      y(5,3)=-1.0
      y(5,6)=-1.0
      y(5,8)=1.0
      y(5,9)=1.0

      y(6,1)=r(1)
      y(6,2)=-r(2)
      y(6,4)=r(4)

      y(7,3)=r(3)
      y(7,6)=-r(6)

      y(8,4)=-r(4)
      y(8,6)=r(6)
      y(8,7)=r(7)
      y(8,8)=r(8)

      y(9,8)=r(8)
      y(9,9)=r(9)
      y(9,10)=r(10)

      y(10,5)=r(5)
      y(10,7)=r(7)

      y(11,2)=r(2)
      y(11,5)=r(5)
      y(11,10)=r(10)

C Generamos un fichero de salida para comprobar que los datos se han
C introducido correctamente:

      open (10,file='2-9-stc.out',status='unknown')
      do i=1,n
         write(10,*) (y(i,j),j=1,n),' | ',b(i)
         write(*,*) (y(i,j),j=1,n),' | ',b(i)
      enddo

      write(10,*)'****Datos correspondientes al ejercicio 2-6-stc.****'
!      write(*,*)    !Imprimimos en pantalla la matriz.

      do i=1,np
        write(10,*) b(i)
        b0(i)=b(i)  !Guardamos los terminos independientes. Que ->
        !correponden con el vector V, y que necesitamos para obtener
        !la resistencia equivalente.
      end do

      close (10)
C Llamamos a las subrutinas para la descomposici¢n y resoluci¢n de la
C matriz
      call ludcmp(y,n,np,indx,d)
      call lubksb(y,n,np,indx,b)
C Con los valores que hemos obtenido ahora para b y el valor anterior
C que hemos guardado podemos calcular la resistencia total del circuito.
C Una explicaci¢n m s detallada de este proceso se encuentra en el
C informe.
      rt=b0(8)/b(8)
C Escribimos el resultado en pantalla
      write(*,*)'La resistencia equivalente es'
      write(*,*) rt
      write(*,*) b0(8)
      write(*,*) b(8)
      write(*,*)'***********ESTE PROGRAMA HA FINALIZADO**********'
      stop
      end

*+++++++++++++++++++++++SUBRUTINAS+++++++++++++++++++++++++++++++++++++*
************************************************************************
************************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
************************************************************************
************************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
************************************************************************
************************************************************************
