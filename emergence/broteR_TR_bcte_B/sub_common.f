
		subroutine AbT(AT,bmT)
		implicit none
      include "common_brote.h"
		double precision AT(2,2), bmT(2), DamC
		integer b1,b2,i,j
		AT=0.d0
		bmT=0.d0
		do b1=0,1
		do b2=0,1
		 do i=1,0,-1
		 do j=1,0,-1
			AT(2-i,2-j)=AT(2-i,2-j)+  x(1,2-b1,2-b2) *( 
     +  				DamC(i,i,j,b1,j,b2,j,b2)-DamC(i,i,1-j,b1,j,b2,j,b2) )
		 enddo
		 bmT(2-i)= bmT(2-i)+  x(1,2-b1,2-b2)*( 
     +	 					 	       x2oG*DamC(i,i,0,b1,1,b2,1,b2)
     +	 				  	  +(1.d0-x2oG)*DamC(i,i,1,b1,0,b2,0,b2)   )						
		 enddo
		enddo
		enddo
		return
		end

c------------------------------------------------------------------------

		double precision function Stt(l1,l2,i1,i2)
		implicit none
      include "common_brote.h"
		integer l1,l2,i1,i2, a1,a2,b1,b2
		double precision DamC
		Stt=0.d0
		do a1=0,1
		do a2=0,1
		do b1=0,1
		do b2=0,1
		 	Stt=Stt+ x(i1,2-a1,2-a2)*x(i2,2-b1,2-b2)
     +										*DamC(l1,l2,a1,b1,a2,b2,a1,b1)
		enddo
		enddo
		enddo
		enddo
		return
		end

		double precision function SttB(l1,l2,i1,i2)
		implicit none
      include "common_brote.h"
		integer l1,l2,i1,i2, a1,a2,b1,b2
		double precision DamC
		SttB=0.d0
		do a1=0,1
		do a2=0,1
		do b1=0,1
		do b2=0,1
		 	SttB=SttB+ x(i1,2-a1,2-a2)*x(i2,2-b1,2-b2)
     +										*DamC(l1,l2,a1,b1,a2,b2,a2,b2)
		enddo
		enddo
		enddo
		enddo
		return
		end

		double precision function Srt(l1,l2,xx1,xx2)
		implicit none
      include "common_brote.h"
		integer l1,l2
		double precision xx1,xx2, Cfun,DamRT
		integer a1,b1
		Srt=0.d0
		do a1=0,1
		do b1=0,1
		 	Srt=Srt+ Cfun(a1,xx1)*Cfun(b1,xx2)*DamRT(l1,l2,a1,b1)
		enddo
		enddo
		return
		end

		double precision function SrtB(l1,l2,xx1,xx2)
		implicit none
      include "common_brote.h"
		integer l1,l2
		double precision xx1,xx2, Cfun,DamRTb
		integer a1,b1
		SrtB=0.d0
		do a1=0,1
		do b1=0,1
		 	SrtB=SrtB+ Cfun(a1,xx1)*Cfun(b1,xx2)*DamRTb(l1,l2,a1,b1)
		enddo
		enddo
		return
		end

		double precision function Str(l1,l2,xx1,xx2)
		implicit none
      include "common_brote.h"
		integer l1,l2
		double precision xx1,xx2, Cfun,DamTR
		integer a1,b1
		Str=0.d0
		do a1=0,1
		do b1=0,1
		 	Str=Str+ Cfun(a1,xx1)*Cfun(b1,xx2)*DamTR(l1,l2,a1,b1)
		enddo
		enddo
		return
		end

		double precision function StrB(l1,l2,xx1,xx2)
		implicit none
      include "common_brote.h"
		integer l1,l2
		double precision xx1,xx2, Cfun,DamTRb
		integer a1,b1
		StrB=0.d0
		do a1=0,1
		do b1=0,1
		 	StrB=StrB+ Cfun(a1,xx1)*Cfun(b1,xx2)*DamTRb(l1,l2,a1,b1)
		enddo
		enddo
		return
		end

		double precision function Sh(xx1,xx2)
		implicit none
      include "common_brote.h"
		double precision xx1,xx2, Cfun,DamH
		integer a1,b1
		Sh=0.d0
		do a1=0,1
		do b1=0,1
		 	Sh=Sh+ Cfun(a1,xx1)*Cfun(b1,xx2)*DamH(a1,b1)
		enddo
		enddo
		return
		end


C==========================================================================

		double precision function DamH(alp,bet)
		implicit none
      include "common_brote.h"
		integer alp,bet
		DamH = (1.d0-eps)*imoral(2-alp,2-bet,2-iaction(2-alp,2-bet))
     +		 + eps*imoral(2-alp,2-bet,2)
		return
		end

		double precision function DamTRb(l1,l2,alp,bet)
		implicit none
      include "common_brote.h"
		integer l1,l2,alp,bet, iCfun
		DamTRb=(1.d0-eps)
     +*iCfun(l1,imoral(2-alp,2-bet,2-iaction(2-alp,2-bet)))
     +*iCfun(l2,1-iaction(2-alp,2-bet))
     +			+eps *iCfun(l1,imoral(2-alp,2-bet,2))
     +					*iCfun(l2,1)
		return
		end

		double precision function DamTR(l1,l2,alp,bet)
		implicit none
      include "common_brote.h"
		integer l1,l2,alp,bet, iCfun
		DamTR=(1.d0-eps)
     +*iCfun(l2,imoral(2-alp,2-bet,2-iaction(2-alp,2-bet)))
     +*iCfun(l1,1-iaction(2-alp,2-bet))
     +			+eps *iCfun(l2,imoral(2-alp,2-bet,2))
     +					*iCfun(l1,1)
		return
		end

		double precision function DamRTb(l1,l2,alp,bet)
		implicit none
      include "common_brote.h"
		integer l1,l2,alp,bet, iCfun
		DamRTb=iCfun(l1,imoral(2-alp,2-bet,2))*iCfun(l2,1)
		return
		end

		double precision function DamRT(l1,l2,alp,bet)
		implicit none
      include "common_brote.h"
		integer l1,l2,alp,bet, iCfun
		DamRT=iCfun(l1,1)*iCfun(l2,imoral(2-alp,2-bet,2))
		return
		end

		double precision function 
     +						DamC(l1,l2,alp1,bet1,alp2,bet2,alpA,betA)
		implicit none
      include "common_brote.h"
		integer l1,l2,alp1,bet1,alp2,bet2,alpA,betA, iCfun
		DamC=(1.d0-eps)
     +*iCfun(l1,imoral(2-alp1,2-bet1,2-iaction(2-alpA,2-betA)))
     +*iCfun(l2,imoral(2-alp2,2-bet2,2-iaction(2-alpA,2-betA)))
     +			+eps *iCfun(l1,imoral(2-alp1,2-bet1,2))
     +					*iCfun(l2,imoral(2-alp2,2-bet2,2))

      return
      end


C***************************************************************************

		subroutine calc_ppB(kp,ap,bp,ppB,den)
		implicit none
      include "common_brote.h"
		double precision kp,ap,bp,ppB, num,den, DamH
		num=ap*(kp*DamH(0,1)+(1.d0-kp)*DamH(0,0))+bp
		den=1.d0 + ap*( -kp*DamH(1,1)-(1.d0-kp)*DamH(1,0)
     +									+kp*DamH(0,1)+(1.d0-kp)*DamH(0,0) )
		if((abs(den).lt.cero).and.(abs(num).gt.cero)) stop 'Error ppB'
		ppB=num/den
		!print*,num,den,DamH(0,0),DamH(0,1),DamH(1,0),DamH(1,1),epsA,epsB
		return
		end

		double precision function pp(ap,bp,PosInd)
		implicit none
      include "common_brote.h"
		integer PosInd
		double precision ap,bp, AA,BB,CC,DamH, solve2
		AA=ap*(DamH(1,1)-DamH(1,0)-DamH(0,1)+DamH(0,0))
		BB = -1.d0+ ap*(DamH(1,0)+DamH(0,1)-2.d0*DamH(0,0))
		CC = ap*DamH(0,0) +bp
		if((abs(AA).lt.cero).and.(abs(BB).lt.cero).and.(PosInd.ne.1))then
			stop 'Error en pp: AA=BB=0'
		else
			pp=solve2(AA,BB,CC)
		endif
		return
		end

c------------ Resuelve 2º grado -----------------------------------------------------

		double precision function solve2(AA,BB,CC)
		implicit none
      include "common_brote.h"
		double precision AA,BB,CC, xx,r,num,den,x1,x2

		AA=(anint(AA*1.d0/cero))*cero
		BB=(anint(BB*1.d0/cero))*cero
		CC=(anint(CC*1.d0/cero))*cero

			if (abs(AA).lt.cero) then
				if (abs(BB).lt.cero) then
					xx=0.5d0  ! sin error en moral estaría degenerado. con error en moral es esto
				else
					xx=-CC/BB
				endif
			elseif (abs(BB).lt.cero) then
				xx=-CC/AA
				if(xx.lt.0.d0) xx=-xx
			elseif (abs(CC).lt.cero) then
				xx=-BB/AA 
				if(AA.gt.0.d0) then
					xx=0.d0
				elseif(xx.lt.0.d0) then
					xx=0.d0
				endif
			else
				r=(BB**2.d0)-(4.d0*AA*CC)
				r=(anint(r*1.d0/cero))*cero
				if (r.lt.0.d0) then
					print*,'Error: raíz de número negativo'
				else
					num=(-BB-(r**0.5d0))
					num=(anint(num*1.d0/cero))*cero
					x1=num/(2.d0*AA)  ! solución más pequeña
					num=(-BB+(r**0.5d0))
					num=(anint(num*1.d0/cero))*cero
					x2=num/(2.d0*AA)

					if((x1.lt.0.d0).or.(x1.gt.1.d0)) then
						if((x2.lt.0.d0).or.(x2.gt.1.d0)) then
							print*,'Error: x_H fuera de rango'
						else
							xx=x2
						endif
					else
						xx=x1 ! la solución más pequeña es la estable
					endif
				endif
			endif
			solve2=(anint(xx*1.d0/cero))*cero

		return
		end


C---------------------------------------------------------------------------

      double precision function Cfun(ip,xf)
      implicit none
      integer ip
      double precision xf
      if (ip.eq.0) then
			Cfun=1.d0-xf
      elseif (ip.eq.1) then
			Cfun=xf
      else
			print*,ip
			stop "Error 01 in Cfun"
      endif
      return
      end

      integer function iCfun(ip,xf)
      implicit none
      integer ip, xf
      if (ip.eq.0) then
			iCfun=1-xf
      elseif (ip.eq.1) then
			iCfun=xf
      else
			print*,ip
			stop "Error 01 in iCfun"
      endif
      return
      end


	subroutine num2actionmoral(inum)
c	GG(0) GB(1) BG(2) BB(3) GGC(4) GGD(5) GBC(6) GBD(7)...	
	implicit none
	integer inum
      include "common_brote.h"
	integer inumt,i,j,k
	iaction=-1
	imoral=-1
	inumt=inum
	do i=1,2
	 do j=1,2
		iaction(i,j)=mod(inumt,2)
		inumt=floor(inumt/2.)
	 enddo
	enddo
	do i=1,2
	 do j=1,2
	  do k=1,2
		imoral(i,j,k)=mod(inumt,2)
		inumt=floor(inumt/2.)
	  enddo
	 enddo
	enddo
	return
	end



c----------------------------------------------------------------------

       SUBROUTINE FINDInv(matrix, inverse, n,nmax, errorflag)
!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n,nmax
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	double precision, INTENT(IN), DIMENSION(nmax,nmax) :: matrix  !Input matrix
	double precision, INTENT(OUT), DIMENSION(nmax,nmax) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	double precision :: m
	double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix
		
	!Augment input matrix with an identity matrix
	DO i = 1, n
	  DO j = 1, 2*n
	  	IF (j <= n ) THEN
	  	  augmatrix(i,j) = matrix(i,j)
	  	ELSE IF ((i+n) == j) THEN
	  	  augmatrix(i,j) = 1
	  	Else
	  	  augmatrix(i,j) = 0
	  	ENDIF
	  END DO
	END DO

	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
	  IF (augmatrix(k,k) == 0) THEN
	  	FLAG = .FALSE.
	  	DO i = k+1, n
	  	  IF (augmatrix(i,k) /= 0) THEN
	  	  	DO j = 1,2*n
	  	  	  augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
	  	  	END DO
	  	  	FLAG = .TRUE.
	  	  	EXIT
	  	  ENDIF
	  	  IF (FLAG .EQV. .FALSE.) THEN
	  	  	PRINT*, "-Matrix is non - invertible"
	  	  	inverse = 0
	  	  	errorflag = -1
	  	  	return
	  	  ENDIF
	  	END DO
	  ENDIF
	  DO j = k+1, n	  	
	  	m = augmatrix(j,k)/augmatrix(k,k)
	  	DO i = k, 2*n
	  	  augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
	  	END DO
	  END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
	  IF (augmatrix(i,i) == 0) THEN
	  	PRINT*, "Matrix is non - invertible"
	  	inverse = 0
	  	errorflag = -1
	  	return
	  ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
	  m = augmatrix(i,i)
	  DO j = i , (2 * n)	  	  
	  	   augmatrix(i,j) = (augmatrix(i,j) / m)
	  END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
	  DO i =1, k
	  m = augmatrix(i,k+1)
	  	DO j = k, (2*n)
	  	  augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
	  	END DO
	  END DO
	END DO	  	  
	
	!store answer
	DO i =1, n
	  DO j = 1, n
	  	inverse(i,j) = augmatrix(i,j+n)
	  END DO
	END DO
	errorflag = 0

      END SUBROUTINE FINDinv


