
		include "sub_common.f"

		program broteR

		implicit none
		include "common_brote.h"
		character*80 outp
		integer i,j,k,stg, iepsB,iepsA, NepsB, errorf,ibb,maymen,
     +								NSOLCmax(Nepsmax),
     +		SENTIDO(Nepsmax,Nbmax,10), encontradouno,NSOLC(Nepsmax,Nbmax)
		double precision epsBini,epsBfin,wepsB,blim, 
     +			calc_xH,calc_x2oG,wepsA,epsAini,epsAfin,
     +			payoff11, payoff21,difpay,difpay1,difpay2,
     +		  payoff12,payoff22,epsA1,epsA2,CURVE(Nepsmax,Nbmax,10)


		print*, 'Number of strategy: '
		read(*,*) stg
		print*, 'b'
		read(*,*) b
		print*,'epsB: initial, final, number: '
		read(*,*) epsBini, epsBfin, NepsB
		print*,'epsA: initial, final, w: '
		read(*,*) epsAini, epsAfin, wepsA
		print*,'Output file'
		read(*,*) outp
						print*,epsAini, epsAfin, wepsA
		wepsB=(epsBfin-epsBini)/(NepsB-1.d0)
		
		Ny=5
       yvect=(/0.01d0,0.3d0,0.5d0,0.7d0,0.99d0,
     &					0.d0,0.d0,0.d0,0.d0,0.d0/)
c		Ny=1
c       yvect=(/0.3d0,0.d0,0.d0,0.d0,0.d0,
c     &					0.d0,0.d0,0.d0,0.d0,0.d0/)
		Nb=4
       bvect=(/1.1d0,1.5d0,2.d0,3.d0,0.d0,0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
c		Nb=1
c      bvect=(/1.5d0,2.d0,3.d0,0.d0,0.d0,0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)

		NSOLC=1
		NSOLCmax=1

		call num2actionmoral(stg)
		print*,'Strategy:  ',stg


      open(90,file=outp,status='UNKNOWN')
	        write(90,'(A1,A)') '#',' epsB vs epsA  y(1st row)'
      !write(90,'(A1,A)') '#',' epsB vs eps   b(1st row)'
		write(90,'(A8,4x,A,$)') '0000',' '
		do i=1,Ny
			!write(90,'(F10.4,$)') bvect(i)
			write(90,'(F10.4,$)') yvect(i)
		enddo
		write(90,'(A,$)')  '      '
		!do i=1,Nb
		!	write(90,'(I4,$)') 0
		!enddo
		do i=1,Ny
			write(90,'(I4,$)') 0
		enddo
		write(90,'(A)') ' '

		SENTIDO=0
		CURVE=99
		!do ibb=1,Nb 
		! b=bvect(ibb)
		! 			print*,'b= ',b
		do ibb=1,Ny 
		 yy=yvect(ibb)
		 			print*,'y= ',yy
		 epsB=epsBini
		 do iepsB=1,NepsB
			encontradouno=0
		 	print*,'epsB= ',epsB
			epsA=epsAini
  43		continue
				call process(payoff11,payoff21,payoff12,payoff22)
				difpay1=payoff11-payoff21
				epsA1=epsA
  41		continue
			epsA=epsA1+wepsA
				call process(payoff11,payoff21,payoff12,payoff22)
				difpay2=payoff11-payoff21
				epsA2=epsA
			if ((difpay1.lt.0.d0).and.(difpay2.lt.0.d0).or.
     +					((difpay1.gt.0.d0).and.(difpay2.gt.0.d0))) then
				if((epsA2.ge.epsAfin).and.(encontradouno.eq.0)) then
					if(difpay1.lt.0.d0) then
						CURVE(iepsB,ibb,1)=1.d0
						SENTIDO(iepsB,ibb,1)=-1.d0
					elseif(difpay1.gt.0.d0) then
						CURVE(iepsB,ibb,1)=0.d0
					endif
c					if ((difpay1.lt.0.d0).and.(difpay2.lt.0.d0)) then
c						if (difpay1.lt.difpay2) then
c							 CURVE(iy,ibb)=1.d0
c							 SENTIDO(iy,ibb)=-1
c						else
c							 CURVE(iy,ibb)=0.d0
c							 SENTIDO(iy,ibb)=1					
c						endif
c					else
c						if (difpay1.gt.difpay2) then
c							 CURVE(iy,ibb)=1.d0
c							 SENTIDO(iy,ibb)=1
c						else
c							 CURVE(iy,ibb)=0.d0
c							 SENTIDO(iy,ibb)=-1					
c						endif
c					endif
						print*, CURVE(iepsB,ibb,1), SENTIDO(iepsB,ibb,1)
				elseif(epsA.le.epsAfin) then
					epsA1=epsA2
					goto 41
				endif
			elseif ((difpay1*difpay2).lt.0.d0) then

  42			continue
								!print*,epsA1,epsA2, difpay1,difpay2
				!epsA=epsA1-difpay1*(epsA1-epsA2)/(difpay1-difpay2)
				epsA=(epsA1+epsA2)/2.d0
				if((abs(epsA2-epsA1).lt.epstol).and.(epsA.le.epsAfin))then
					if(encontradouno.ge.1) then
							NSOLC(iepsB,ibb)=NSOLC(iepsB,ibb)+1
							NSOLCmax(iepsB)=NSOLC(iepsB,ibb)
					endif
					CURVE(iepsB,ibb,NSOLC(iepsB,ibb))=epsA
					if(difpay1.lt.difpay2) SENTIDO(iepsB,ibb,NSOLC(iepsB,ibb))=-1
					if(difpay2.lt.difpay1) SENTIDO(iepsB,ibb,NSOLC(iepsB,ibb))=1
								print*, CURVE(iepsB,ibb,NSOLC(iepsB,ibb)),
     +							SENTIDO(iepsB,ibb,NSOLC(iepsB,ibb)),NSOLC(iepsB,ibb)
					encontradouno=1
					if(epsA2.lt.epsAfin) then
						epsA=epsA2
						goto 43
					endif
				elseif(epsA.le.epsAfin) then
					call process(payoff11,payoff21,payoff12,payoff22)
					difpay=payoff11-payoff21
										
							!print*,epsA,b,difpay

					if(difpay.eq.0.d0) difpay=cero
					if((difpay*difpay1).gt.0.d0) then
						difpay1=difpay
						epsA1=epsA
						goto 42
					elseif((difpay*difpay2).gt.0.d0) then
						difpay2=difpay
						epsA2=epsA
						goto 42
					endif

				endif

			endif

			epsB=epsB+wepsB
		 enddo !epsB


		enddo !b


		epsB=epsBini
		do iepsB=1,NepsB
		 do i=1,NSOLCmax(iepsB)
			write(90,'(F8.4,4x,A,$)') epsB,' '
			do ibb=1,Ny
				write(90,'(F10.4,$)') CURVE(iepsB,ibb,i)
			enddo
			write(90,'(A,$)') '      '
			do ibb=1,Ny
				write(90,'(I4,$)') SENTIDO(iepsB,ibb,i)
			enddo
			write(90,'(A)') ' '
		 enddo
		 epsB=epsB+wepsB
		enddo


		close(90)

		stop
		end




		subroutine process(payoff11,payoff21,payoff12,payoff22)		
		implicit none
		include "common_brote.h"
		double precision payoff11,payoff21,payoff12,payoff22
		integer errorf
		double precision den,pp

					!epsB=epsA

				eps=epsA
				x1AGo=pp(yy,1.d0-yy,0)
				eps=epsB	
				x1BGo=pp(1.d0-yy,yy,0)			
				call IC
				call solve_x1
				x1BoG=x(3,1,1)+x(3,2,1)
				x1AoG=x(1,1,1)+x(1,2,1)
				eps=epsA
				call calc_ppB(yy*x1AoG+(1.d0-yy)*x1BoG,1.d0,0.d0,x2oG,den)
				if(abs(den).lt.cero) stop 'den=0 en ppal'
				call solve_x2(errorf)
				eps=epsA
				x2Go=x(2,1,1)+x(2,1,2)
				call calcpay(payoff11,payoff21,payoff12,payoff22)
					!print*, epsA,payoff11-payoff21, payoff12-payoff22
					!print*,payoff11,payoff21!,payoff12,payoff22
			!print*,x(1,1,1),x(1,1,2) 
			!print*,x(1,2,1),x(1,2,2)
			!print*,x(1,1,1)+x(1,2,2)+x(1,1,2)+x(1,2,1)
			!print*,x(3,1,1),x(3,1,2) 
			!print*,x(3,2,1),x(3,2,2)
			!print*,x(3,1,1)+x(3,2,2)+x(3,1,2)+x(3,2,1)
			!print*
		return
		end


c----------- Solves ---------------------------------------------------------------

		subroutine solve_x1

      implicit none
      include "common_brote.h"
		integer a1,a2, cont, is
		double precision difx2(3),xGGold(3),xGBold(3), suma(3),
     +		xBGold(3),xBBold(3), x11new,x22new,Stt,Srt,Sh,x12new,x21new

		cont=0
		difx2=1.d0
		x1AoG=x(1,1,1)+x(1,2,1)
		x1BoG=x(3,1,1)+x(3,2,1)

						!print*
						!print*,x(1,1,1),x(1,1,2)
						!print*,x(1,2,1),x(1,2,2)
						!print*,'---',x1AGo
						!print*,x(3,1,1),x(3,1,2)
						!print*,x(3,2,1),x(3,2,2)
						!print*,'---',x1BoG
						!print*

		do while ((difx2(1)**0.5d0.gt.xtol).or.(difx2(3)**0.5d0.gt.xtol))
		 cont=cont+1

		 do is=1,3,2
		  xGGold(is)=x(is,1,1)
		  xGBold(is)=x(is,1,2)
		  xBGold(is)=x(is,2,1)
		  xBBold(is)=x(is,2,2)
		 enddo

		  eps=epsA
		  x11new= dt*(yy*Stt(1,1,1,1)+(1.d0-yy)*Srt(1,1,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(1,1,1)
		  x22new= dt*(yy*Stt(0,0,1,1)+(1.d0-yy)*Srt(0,0,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(1,2,2)
		  x12new= dt*(yy*Stt(1,0,1,1)+(1.d0-yy)*Srt(1,0,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(1,1,2)
		  x21new= dt*(yy*Stt(0,1,1,1)+(1.d0-yy)*Srt(0,1,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(1,2,1)
		  x(1,1,1)=x11new
		  x(1,2,2)=x22new
		  x(1,1,2)=x12new
		  x(1,2,1)=x21new
		  !x(1,1,2)=x1AGo-x11new
		  !x(1,2,1)=1.d0-x1AGo-x22new


		  eps=epsB
		  x11new= dt*((1.d0-yy)*Stt(1,1,3,3)+yy*Srt(1,1,x1BoG,x1AoG))
     +												 +(1.d0-dt)*x(3,1,1)
		  x22new= dt*((1.d0-yy)*Stt(0,0,3,3)+yy*Srt(0,0,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(3,2,2)
		  x12new= dt*((1.d0-yy)*Stt(1,0,3,3)+yy*Srt(1,0,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(3,1,2)
		  x21new= dt*((1.d0-yy)*Stt(0,1,3,3)+yy*Srt(0,1,x1AoG,x1BoG))
     +												 +(1.d0-dt)*x(3,2,1)
		  x(3,1,1)=x11new
		  x(3,2,2)=x22new
		  x(3,1,2)=x12new
		  x(3,2,1)=x21new
		  !x(3,1,2)=x1BGo-x11new
		  !x(3,2,1)=1.d0-x1BGo-x22new

		  do is=1,3,2
			do a1=1,2
			do a2=1,2
				if(x(is,a1,a2).lt.0.d0) x(is,a1,a2)=0.d0
			enddo
			enddo
			suma(is)=x(is,1,1)+x(is,2,1)+x(is,1,2)+x(is,2,2)
			if(suma(is).ne.1.d0) then
			do a1=1,2
			do a2=1,2
				x(is,a1,a2)=x(is,a1,a2)/suma(is)
			enddo
			enddo	
			endif
		  enddo

		  x1AoG=x(1,1,1)+x(1,2,1)
		  x1BoG=x(3,1,1)+x(3,2,1)

			!print*,x(1,1,1),x(1,1,2) 
			!print*,x(1,2,1),x(1,2,2)
			!print*,x(1,1,1)+x(1,2,2)+x(1,1,2)+x(1,2,1)
			!print*,x(3,1,1),x(3,1,2) 
			!print*,x(3,2,1),x(3,2,2)
			!print*,x(3,1,1)+x(3,2,2)+x(3,1,2)+x(3,2,1)
			!print*
			!if(cont.eq.10) stop

		 do is=1,3,2
		  difx2(is)=(x(is,1,1)-XGGold(is))**2.d0
     +		+(x(is,1,2)-XGBold(is))**2.d0
     +      +(x(is,2,1)-XBGold(is))**2.d0+(x(is,2,2)-XBBold(is))**2.d0
		 enddo
		enddo

		do is=1,3,2
		do a1=1,2
		do a2=1,2
			x(is,a1,a2)=(anint(x(is,a1,a2)*1.d0/cero))*cero
		enddo
		enddo
		enddo
				!stop
		return 
		end




		subroutine solve_x2(errorf)

      implicit none
      include "common_brote.h"
		integer i,j,b1,b2,errorf
		double precision  AT(2,2),bmT(2),AR(2,2),bmR(2),A(2,2),bm(2),
     +						Ainv(2,2),xsol(2), DamRTb, Str

		call AbT(AT,bmT)

		AR=0.d0
		bmR=0.d0
		do i=1,0,-1
		 	bmR(2-i)= Str(i,i,x2oG,x1BoG)
		enddo


		do i=1,2
		do j=1,2
			A(i,j)=yy*AT(i,j)
		enddo
		A(i,i)=A(i,i)-1.d0
		bm(i)=yy*bmT(i)+(1.d0-yy)*bmR(i)
		enddo
		bm=-bm

		call FINDInv(A,Ainv,2,2,errorf)
		if (errorf.ne.0) then
			goto 51
		endif
		xsol=matmul(Ainv,bm)
		x(2,1,1)=xsol(1)
		x(2,2,2)=xsol(2)
		x(2,2,1)=x2oG-x(2,1,1)
		x(2,1,2)=1.d0-x2oG-x(2,2,2)

 51	continue
		do b1=1,2
		do b2=1,2
			x(2,b1,b2)=(anint(x(2,b1,b2)*1.d0/cero))*cero
			if(x(2,b1,b2).lt.(-1.d0*cero)) then
				print*,b1,b2,x(2,b1,b2)
				print*, 'menor que 0 ***********************'
				!stop
			endif
		enddo
		enddo
		return 
		end



		subroutine IC
		implicit none
      include "common_brote.h"
		integer i
		double precision xHA,xHB, pp

		eps=epsA
		xHA=pp(1.d0,0.d0,1)
			!print*,'xHA= ',xHA,epsA
		eps=epsB
		xHB=pp(1.d0,0.d0,1)
			!print*,'xHB= ',xHB,epsA
			!stop
		do i=3,3
			x(i,1,1)=xHB
			x(i,2,2)=1.d0-xHB
			x(i,1,2)=0.d0
			x(i,2,1)=0.d0
		enddo
		do i=1,2
			x(i,1,1)=xHA
			x(i,2,2)=1.d0-xHA
			x(i,1,2)=0.d0
			x(i,2,1)=0.d0
		enddo
		return
		end



c------------ payoff -----------------------------------------------------------------


		subroutine calcpay(payoff11,payoff21,payoff12,payoff22)

		implicit none
      include "common_brote.h"
		double precision  payoff11, payoff21,
     +						payoff12,payoff22
		double precision p11,p12,p21,p22, Cfun,num,den
		integer ia, ib

		p11=0.d0
		p12=0.d0
		p21=0.d0
		p22=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
			p11=p11+ Cfun(ia,x1AGo)*Cfun(ib,x1AGo)*iaction(2-ia,2-ib)  		! 1 le da a 1
			p21=p21+ Cfun(ia,x2oG)*( yy*Cfun(ib,x1AoG) 								! 2 le da a 1
     +              +(1.d0-yy)*Cfun(ib,x1BoG)	)*iaction(2-ia,2-ib)				 
			p12=p12+ Cfun(ia,x1AGo)*Cfun(ib,x2Go) *iaction(2-ia,2-ib)			! 1 le da a 2 
			p22=p22+ Cfun(ia,x2oG)*Cfun(ib,x2oG)*iaction(2-ia,2-ib)           ! 2 le da a 2  comprobar
		enddo
		enddo
		p11=yy*p11
		p12=yy*p12

			payoff11=(b-c)*p11  *(1.d0-epsA)    ! de 1 (con 1)
			payoff21=(b*p12-c*p21)  *(1.d0-epsA)    ! de 2 (con 1)
			payoff12=(b*p21-c*p12)  *(1.d0-epsA)    ! de 1 (con 2)
			payoff22=(b-c)*p22   *(1.d0-epsA)	 ! de 2 (con 2)

		return
		end
