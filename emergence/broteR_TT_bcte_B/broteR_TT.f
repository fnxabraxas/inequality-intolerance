
		include "sub_common.f"

		program broteR

		implicit none
		include "common_brote.h"
		character*80 outp
		integer i,j,k,stg, iepsB,iepsA, NepsB, errorf,ibb,maymen,
     +								NSOLCmax(Nepsmax),
     +		SENTIDO(Nepsmax,Nbmax,11), encontradouno,NSOLC(Nepsmax,Nbmax)
		double precision epsBini,epsBfin,wepsB,blim, 
     +			calc_xH,calc_x2oG,wepsA,epsAini,epsAfin,
     +			payoff11, payoff21,difpay,difpay1,difpay2,
     +		  payoff12,payoff22,epsA1,epsA2,CURVE(Nepsmax,Nbmax,11)


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
		
		Ny=11
       yvect=(/0.01d0,0.1d0,0.2d0,0.3d0,0.4d0,
     &				0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,0.99d0/)
c		Ny=1
c       yvect=(/0.3d0,0.d0,0.d0,0.d0,0.d0,
c     &					0.d0,0.d0,0.d0,0.d0,0.d0/)
c		Nb=4
c       bvect=(/1.1d0,1.5d0,2.d0,3.d0,0.d0,0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
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
	
					!epsA=0.3
					!epsB=0.7
					!epsB=epsA
				call solve_x1AB(x1AGo,x1BGo)
					!print*, 'eps= ', epsA,epsB
					!print*,x1AGo,x1BGo
				call IC
					!print*,'IC: ',x
				eps=epsA
				call solve_x1A
				x1AoG=x(1,1,1)+x(1,2,1)
					!print*,'x1AoG= ',x1AoG,x(1,1,1),x(1,2,1)
				call calc_ppB(x1AoG,yy,1.d0-yy,x2oG,den)
				if(abs(den).lt.cero) stop 'den=0 en ppal'
				call solve_x2(errorf)
				x2Go=x(2,1,1)+x(2,1,2)
					!print*,'x2Go= ',x2Go
				call calcpay(payoff11,payoff21,payoff12,payoff22)
			!print*,'Pays: ',payoff11,payoff21 !,payoff12,payoff22
					!stop
		return
		end


c----------- Solves ---------------------------------------------------------------


		subroutine solve_x1AB(xA,xB)

      implicit none
      include "common_brote.h"
		double precision xA,xB
		integer a1,a2, cont
		double precision difx2,xAold,xBold,xAch,xBch,Cfun,DamH,pp

						!print*,x(1,1,1),x(1,1,2), Stt(1,1,1,1)
						!print*,x(1,2,1),x(1,2,2), Stt(0,0,1,1)
						!print*,'---',x1AGo

				eps=epsA
				xA=pp(1.d0,0.d0,1)
				eps=epsB
				xB=pp(1.d0,0.d0,1)
					!print*,'xA,xB= ',xA,xB

		cont=0
		difx2=1.d0
		do while (difx2**0.5d0.gt.xtol)
		  cont=cont+1

		  xAold=xA
		  xBold=xB

		  eps=epsA
		  xAch=0.d0
		  do a1=0,1
		  do a2=0,1
		 	xAch=xAch+yy*Cfun(a1,xA)*Cfun(a2,xA)*DamH(a1,a2)
     +			  +(1.d0-yy)*Cfun(a1,xA)*Cfun(a2,xB)*DamH(a1,a2)
		  enddo
		  enddo
		  xA=(1.d0-dt)*xA+dt*xAch
		  eps=epsB
		  xBch=0.d0
		  do a1=0,1
		  do a2=0,1
		 	xBch=xBch+yy*Cfun(a1,xB)*Cfun(a2,xA)*DamH(a1,a2)
     +			  +(1.d0-yy)*Cfun(a1,xB)*Cfun(a2,xB)*DamH(a1,a2)
		  enddo
		  enddo		 
		  xB=(1.d0-dt)*xB+dt*xBch

		  difx2=(xA-xAold)**2.d0+(xB-xBold)**2.d0
		enddo

		return 
		end




		subroutine solve_x1A

      implicit none
      include "common_brote.h"
		integer a1,a2, cont, is
		double precision difx2, xGGold, xGBold,xBGold, xBBold,suma(3),
     +							 x11new,x22new,x12new,x21new, Stt,StrB

						!print*,x(1,1,1),x(1,1,2), Stt(1,1,1,1)
						!print*,x(1,2,1),x(1,2,2), Stt(0,0,1,1)
						!print*,'---',x1AGo

		cont=0
		difx2=1.d0
		do while (difx2**0.5d0.gt.xtol)
		  cont=cont+1

		  xGGold=x(1,1,1)
		  xGBold=x(1,1,2)
		  xBGold=x(1,2,1)
		  xBBold=x(1,2,2)
		  x11new= dt*(yy*Stt(1,1,1,1)+(1.d0-yy)*StrB(1,1,x1AGo,x1BGo))
     +												 +(1.d0-dt)*x(1,1,1)
		  x22new= dt*(yy*Stt(0,0,1,1)+(1.d0-yy)*StrB(0,0,x1AGo,x1BGo))
     +												 +(1.d0-dt)*x(1,2,2)
		  x12new= dt*(yy*Stt(1,0,1,1)+(1.d0-yy)*StrB(1,0,x1AGo,x1BGo))
     +												 +(1.d0-dt)*x(1,1,2)
		  x21new= dt*(yy*Stt(0,1,1,1)+(1.d0-yy)*StrB(0,1,x1AGo,x1BGo))
     +												 +(1.d0-dt)*x(1,2,1)
		  x(1,1,1)=x11new
		  x(1,2,2)=x22new
		  x(1,1,2)=x12new
		  x(1,2,1)=x21new

		  do is=1,1!3,2
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

			!print*,x(1,1,1),x(1,1,2), Stt(1,1,1,1), StrB(1,1,x1AGo,x1AGo)
			!print*,x(1,2,1),x(1,2,2), Stt(0,0,1,1), StrB(0,0,x1AGo,x1AGo)
			!print*
			!if(cont.eq.10) stop

		  difx2=(x(1,1,1)-XGGold)**2.d0+(x(1,1,2)-XGBold)**2.d0
     +		+(x(1,2,1)-XBGold)**2.d0+(x(1,2,2)-XBBold)**2.d0
		enddo

		do a1=1,2
		do a2=1,2
			x(1,a1,a2)=(anint(x(1,a1,a2)*1.d0/cero))*cero
		enddo
		enddo

		return 
		end




		subroutine solve_x2(errorf)

      implicit none
      include "common_brote.h"
		integer i,j,b1,b2,errorf
		double precision  AT(2,2),bmT(2),AR(2,2),bmR(2),A(2,2),bm(2),
     +						Ainv(2,2),xsol(2), DamC,DamRTb, Cfun

		call AbT(AT,bmT)

		AR=0.d0
		bmR=0.d0
		do b1=0,1
		 do i=1,0,-1
		 do j=1,0,-1
			AR(2-i,2-j)=AR(2-i,2-j)+ Cfun(b1,x1AGo) *(
     +					DamRTb(i,i,j,b1)-DamRTb(i,i,1-j,b1)	)
		 enddo
		 bmR(2-i)= bmR(2-i)+ Cfun(b1,x1AGo) *( 
     +	 					 	       x2oG*DamRTb(i,i,0,b1)
     +	 				  	  +(1.d0-x2oG)*DamRTb(i,i,1,b1)   )	
		 enddo
		enddo


		do i=1,2
		do j=1,2
			A(i,j)=yy*AT(i,j)+(1.d0-yy)*AR(i,j)
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
		do i=1,3
			x(i,1,1)=x1AGo
			x(i,2,2)=1.d0-x1AGo
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
		double precision p12,p21,p22, Cfun, p11A, p1A1
		integer ia, ib

		p11A=0.d0
		p1A1=0.d0
		p12=0.d0
		p21=0.d0
		p22=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
			p11A=p11A+ Cfun(ib,x1AGo)*( yy*Cfun(ia,x1AGo)*(1.d0-epsA) 
     &	 + (1.d0-yy)*Cfun(ia,x1BGo)*(1.d0-epsB) )*iaction(2-ia,2-ib)  		! 1 le da a 1A
			p1A1=p1A1+ Cfun(ia,x1AGo)*( yy*Cfun(ib,x1AGo)*(1.d0-epsA) 
     &	 + (1.d0-yy)*Cfun(ib,x1BGo)*(1.d0-epsA) )*iaction(2-ia,2-ib)  		! 1A le da a 1
			p22=p22+ Cfun(ia,x2oG)*Cfun(ib,x2oG)*iaction(2-ia,2-ib) 	! 2 le da a 2
			p21=p21+ Cfun(ia,x2oG)*Cfun(ib,x1AoG)*iaction(2-ia,2-ib)	! 2 le da a 1
			p12=p12+ Cfun(ia,x1AGo)*Cfun(ib,x2Go)*iaction(2-ia,2-ib)   	! 1 le da a 2
		enddo
		enddo
		p21=yy*p21 *(1.d0-epsA) 
		p22=p22 *(1.d0-epsA)
		p12=p12*(yy*(1.d0-epsA)+(1-yy)*(1.d0-epsB))

				!print*,p22,p21,p12,p1A1,p11A

			payoff11=b*p11A -c*p1A1   ! de 1 (con 1)
			payoff21=(b*p12-c*p21)   ! de 2 (con 1)
			payoff12=(b*p21-c*p12)   ! de 1 (con 2)
			payoff22=((b-c)*p22  ) 	 ! de 2 (con 2)

				!print*, payoff11, payoff21, epsA, epsB


		return
		end



		subroutine limb(blim,maymen)
c		maymen:  +1: b>   -1: b<	para que haya invasión	

		implicit none
      include "common_brote.h"
		double precision  blim
		integer maymen
		double precision p11,p12,p21,p22, Cfun,num,den
		integer ia, ib

		p11=0.d0
		p12=0.d0
		p21=0.d0
		p22=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
			p11=p11+ Cfun(ia,x1AGo)*Cfun(ib,x1AGo)*iaction(2-ia,2-ib)  		! 1 le da a 1
			p22=p22+ Cfun(ia,x2oG)*Cfun(ib,x2oG)*iaction(2-ia,2-ib) 	! 2 le da a 2
			p21=p21+ Cfun(ia,x2oG)*Cfun(ib,x1AoG)*iaction(2-ia,2-ib)	! 2 le da a 1
			p12=p12+ Cfun(ia,x1AGo)*Cfun(ib,x2Go)*iaction(2-ia,2-ib)   	! 1 le da a 2
		enddo
		enddo
		p21=yy*p21

		num=p11-p21
		den=p11-p12
		blim=num/den
			print*,num,den
		if(den.gt.1.d-12) then
			maymen=-1
		elseif(den.lt.-1.d-12) then
			maymen=1
		else
			maymen=0
			blim=1.d0
		endif

		return
		end






