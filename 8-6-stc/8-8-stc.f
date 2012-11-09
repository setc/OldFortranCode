************************************************************************
*     8-6a-stc                                                         *
*     Calculo de la transformada de Fourier R pida para la concentraci-*
*     ¢n de CO2 en Hawai 1959-2002, teniendo en cuenta tanto el        *
*     aumento sinuoidal como el lineal.                                *
*     Sebasti n Torrente Carrillo                                      *
************************************************************************

!En este caso dimensionaremos por defecto, tomando 2**9
*      parameter(nmax=2**10)  ! Dimensiono por exceso
      parameter(nmax=2**9)   ! Dimensiono por defecto
      dimension fr(0:nmax-1),fi(0:nmax-1),data(1:2*nmax),
     &          gr(0:nmax-1),gi(0:nmax-1)
      write(*,*) 'nmax=',nmax

* Creamos una variable contador y la inicializamos a cero:
      k=0

* Declaramos todos los valores de la funci¢n a cero.
      do j=1,nmax-1
        fr(j)=0.0   ! Parte real
        fi(j)=0.0   ! Parte imaginaria
      end do

      open(10,file='co2hawaic.dat',blank='null')   ! fichero de lectura
      open(11,file='8-6a-stc.dat',status='unknown')    ! fichero de escritura
      read(10,*)     ! salto la lectura de la primera fila (que es encabezamiento)
* leo los 12 meses de cada an¤o, desde 1959 hasta 2002

      do i=1,nmax            ! esto es cuando dimensiono por defecto (pierdo datos)
*      do i=1,(2002-1959)*12  ! esto es cuando dimensiono por exceso (introduzco datos nulos)

* Mes (m) y la concentraci¢n de co2:
         read (10,*) m,co2c
* Y escribo cada par de variables en el '8-6a-stc.dat'
         write (11,*) m,co2c
* defino un contador que sirve de indice a f_k
         k=k+1
         fr(k)=co2c   ! parte real de la senyal muestreada
         fi(k)=0.0   ! parte imaginaria de la senyal muestreada
      end do
      write(*,*) 'numero de meses=',(2002-1959)*12
      close (10)
      close (11)

* Ponemos como numero m ximo de meses los leidos, si dimensionamos por
* defecto:
      tmax=nmax
* Y cuando dimensiono por exceso (datos de exceso=0.0)
*      tmax=(2004-1749)*12.0

      dt=1.0      ! El paso es mes a mes.
      dv=1.0/tmax

* Insertamos en una secuencia de data nuestros datos.
      do i=1,2*nmax,2
        ind=(i-1)/2
        data(i)=fr(ind)
        data(i+1)=fi(ind)
      end do

* isign=1 reemplaza data(1:2*nmax) por su transf. de Fourier discreta
* isign=-1 reemplaza data(1:2*nmax) por nn veces su transf. de Fourier
* inversa discreta
      isign=1
      call four1(data,nmax,isign)

      do i=1,2*nmax,2
        ind=(i-1)/2
        gr(ind)=data(i)
        gi(ind)=data(i+1)
      end do

* Y ahora calculamos la transformada de Fourier:

* Escribimos los datos que generamos
*      open(30,file='8-6a-stc.dat')   ! cuando dimensiono por exceso
      open(30,file='8-6-stc.dat')   ! cuando dimensiono por defecto
      do j=1,nmax-1
        t=dt*j
        write(30,*) t,fr(j),fi(j)
      end do
      close(30)

* Calculamos el espectro de potencial de las frecuencias positivas, sin
* superar la frecuencia de Nyquist.
*      open(40,file='8-6g9-stc.dat') ! cuando dimensiono por exceso
       open(40,file='8-6g9-stc.dat') ! cuando dimensiono por defecto
      do j=1,nmax/2-1
        v=dv*j
        gg=sqrt(gr(j)**2+gi(j)**2)
        write(40,*) v,gg
      end do
      close(40)
      write(*,*) 'frec. de Nyquist =',1.0/(2.0*dt)
      write(*,*) 'programa finalizado'
      stop
      end

************************************************************
      SUBROUTINE four1(data,nn,isign)
* subrutina four1.for del Numerical Recipes (2a ed.) p.501
************************************************************
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
