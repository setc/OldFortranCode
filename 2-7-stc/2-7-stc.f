************************************************************************
*     2-7-stc                                                          *
*     Calculo de una resistencia equivalente, con una resistencia      *
*     variable dentro del sistema                                      *
*     Sebasti n Torrente Carrillo                                      *
*     31-3-2005                                                        *
************************************************************************

C Como se puede apreciar este programa es muy similar al anterior
C (incluso muchas partes de este programa son recortes del anterior)
C esto se debe a la similaridad entre ambos problemas. La £nica
C diferencia es que aqu¡ una de la resistencias es variable.

C En este programa calcularemos la resistencia equivalente de una red
C de resistencias dadas empleando las reglas de Kirchoff y los metodos
C de resoluci¢n de ecuaciones lineales.

C Declaramos variables y procedemos a dar valores a los elementos de la
C matriz. Existe una diferencia con respecto al anterior ejercicio
C y es que una de las matrices tiene un valor variable

C Declaramos la dimensi¢n de la matriz juntoa otras variables.
      integer i,j,np,n,k  ! i y j son indices para los bucles de las
      !matrices, k es para el bucle que incluye a las subrutinas.
      parameter (n=13,np=13)
      dimension r(n-1),y(n,n),b(n),b0(n),indx(np)
      real v!,dk ! correponde al paso dentro del bucle en el que se
      !encuentran las subrutinas.

C Damos valores a v y las resistencias. La £nica excepci¢n es la
C resistencia variable, cuyos valores ser n declarados dentro del
C bucle.

      v=10.0 !Valor arbitrario.
      res=1.0


C Inicializamos todos los valores de la matriz a cero
      do i=1,n
        b(i)=0.0
        r(i)=0.0
        do j=1,n
           y(i,j)=0.0
        enddo
      enddo

C Y ahora asignamos valores:
      b(13)=v

      r(1)=res !Tomamos 1.0 como valor inicial, dentro del bucle este
      !valor se ir  incrementando hasta 100.
      r(2)=res
      r(3)=res
      r(4)=res
      r(5)=res
      r(6)=res
      r(7)=res
      r(8)=res
      r(9)=res
      r(10)=res
      r(11)=res
      r(12)=res
      
C Definimos los valores de la matriz distintos de cero:
      y(1,1)=1.0
      y(1,4)=1.0
      y(1,7)=1.0
      y(1,13)=-1.0
      
      y(2,1)=1.0
      y(2,2)=-1.0
      y(2,8)=-1.0
      
      y(3,3)=-1.0
      y(3,4)=1.0
      y(3,8)=-1.0
      
      y(4,2)=1.0
      y(4,3)=1.0
      y(4,6)=-1.0

      y(5,7)=-1.0
      y(5,9)=-1.0
      y(5,10)=1.0
      
      y(6,5)=-1.0
      y(6,10)=-1.0
      y(6,11)=1.0
      
      y(7,8)=-1.0
      y(7,9)=1.0
      y(7,12)=1.0
      
C     y(8,1)= Se declarar  dentro del bucle, ya que tiene valor variable
      y(8,2)=r(2)
      y(8,3)=-r(3)
      y(8,4)=-r(4)
      
      y(9,4)=-r(4)
      y(9,7)=r(7)
      y(9,9)=-r(9)
      
C     y(10,1)
      y(10,5)=r(5)
      y(10,7)=r(7)
      y(10,10)=-r(10)
      
      y(11,2)=-r(2)
      y(11,5)=r(5)
      y(11,6)=-r(6)
      y(11,11)=r(11)
      
      y(12,9)=r(9)
      y(12,10)=r(10)
      y(12,11)=r(11)
      y(12,12)=-r(12)
      
C      y(13,1)
      y(13,5)=r(5)
      y(13,11)=r(11)

C Generamos un fichero de salida para comprobar que los datos se han
C introducido correctamente, eso s¡, los datos correspondientes a la
C resistencia variable no figuran, por lo que en esos lugares solo
C tendremos el valor 0, de todos modos podemos revisar si existen
C errores.

      open (10,file='2-7a-stc.out',status='unknown')
      do i=1,n
         write(10,*) (y(i,j),j=1,n),' | ',b(i)
         write(*,*) (y(i,j),j=1,n),' | ',b(i)
      enddo

      write(10,*)'****Datos correspondientes al ejercicio 2-7-stc.****'
      close (10)


      do i=1,np
!        write(10,*) b(i)
        b0(i)=b(i)  !Guardamos los terminos independientes. Que ->
        !correponden con el vector V, y que necesitamos para obtener
        !la resistencia equivalente.
      end do


C Llamamos a las subrutinas para la descomposici¢n y resoluci¢n de la
C matriz, en este caso las llamaremos dentro de un bucle, para
C obtener los distintos valores de la resistencia equivalente. Es ahora
C (y dentro del bucle) cuando declaramos los valores de la matriz
C que contienen a la resistencia r(1)
C Antes abrimos el documento de salida que contendr  los datos.
      open (20,file='2-7b-stc.out',status='unknown')
      write(20,*)'#r(1)  r_equiv'
      do k=1,100
        call ludcmp(y,n,np,indx,d)
        call lubksb(y,n,np,indx,b)
        r(1)=r(1)+1.0 !No hemos puesto la obtenci¢n del paso porque ya
                      !hemos decidido de antemano que este sea de uno
                      !en uno.
C Con los valores que hemos obtenido ahora para b y el valor anterior
C que hemos guardado podemos calcular la resistencia total del circuito.
C Una explicaci¢n m s detallada de este proceso se encuentra en el
C informe.
        rt=b0(13)/b(13)
C Y escribimos el resultado en pantalla
        write(*,*) rt
C Y en el documento para poder representarlo graficamente con GNUplot.
        write(20,*)r(1),'   ',rt
      enddo
      close(20)
C Y damos por finalizado el programa.
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
