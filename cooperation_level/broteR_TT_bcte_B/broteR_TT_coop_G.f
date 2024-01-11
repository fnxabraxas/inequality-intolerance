c Calculates cooperation levels before invasion
c Each row ia different epsA. each column is a different epsB
c Results are neither multiplied by (1-epsA) nor y or (1-y): 
c   coopR_TOT = ( y*coopRwA + (1-y)*coopRwB ) *(1-epsA)
c   coopM_TOT = ( y*coopMwA ) *(1-epsA) 

		include "sub_common.f"

		program broteR

		implicit none
		include "common_brote.h"
		character*80 outp,
     +		  outfileARtoA, outfileARtoB, outfileAMtoA,
     + 		outfileBRtoA, outfileBRtoB
		integer i,j,k,stg, NG, iG
		double precision Gini,Gfin, epsT, GG, wG,
     +		 coopARwA,coopARwB,coopBRwA,coopBRwB, coopMwA,
     +			GiniA, GiniB, GfinA, GfinB


		print*, 'Number of strategy: '
		read(*,*) stg
		print*, 'Proportion of A (y): '
		read(*,*) yy
		print*, 'epsT: '
		read(*,*) epsT
		print*,'number of G'
		read(*,*) NG
	print*, 'Output files (5): '
	read(*,*) outfileARtoA, outfileARtoB, outfileBRtoA,
     +			 outfileBRtoB, outfileAMtoA  
						!print*,epsAini, epsAfin, wepsA
	GfinB = (1.d0-yy)
	GfinA = yy*epsT/(1.d0-epsT)
	if (GfinB<GfinA) then
		Gfin=GfinB
	else
		Gfin=GfinA
	endif
	GiniB = -(1.d0-yy)*epsT/(1.d0-epsT)
	GiniA = -yy
	if (GiniB>GiniA) then
		Gini=GiniB
	else
		Gini=GiniA
	endif				
	wG=(Gfin-Gini)/(NG-1.d0)

		call num2actionmoral(stg)
		print*,'Strategy, y:  ',stg,yy


      
      open(91,file=outfileARtoA,status='UNKNOWN')
      open(92,file=outfileARtoB,status='UNKNOWN')
      open(93,file=outfileBRtoA,status='UNKNOWN')
      open(94,file=outfileBRtoB,status='UNKNOWN')
      open(95,file=outfileAMtoA,status='UNKNOWN')
      
      GG=Gini
      do iG=1,NG

	epsA = epsT - (1.d0-epsT)*GG/yy
	epsB = epsT + (1.d0-epsT)*GG/(1.d0-yy)
	!print*,Gini,gfin,epsT,yy,GG,epsA,epsB
      	call process(coopARwA,coopARwB,coopBRwA,coopBRwB, coopMwA)          
        write(91,'(2F10.4,4x,2F10.4)') GG, coopARwA, epsA, epsB
        write(92,'(2F10.4,4x,2F10.4)') GG, coopARwB, epsA, epsB
        write(93,'(2F10.4,4x,2F10.4)') GG, coopBRwA, epsA, epsB
        write(94,'(2F10.4,4x,2F10.4)') GG, coopBRwB, epsA, epsB
        write(95,'(2F10.4,4x,2F10.4)') GG, coopMwA, epsA, epsB
        GG=GG+wG

      enddo !iG
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      
      return
      end


	subroutine process(coopARwA,coopARwB,coopBRwA,coopBRwB, coopMwA) !(payoff11,payoff21,payoff12,payoff22)		
	implicit none
	include "common_brote.h"
	double precision coopARwA,coopARwB,coopBRwA,coopBRwB, coopMwA !payoff11,payoff21,payoff12,payoff22
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
				!call calcpay(payoff11,payoff21,payoff12,payoff22)
			!print*,'Pays: ',payoff11,payoff21 !,payoff12,payoff22
					!stop
	call calccoop(coopARwA,coopARwB,coopBRwA,coopBRwB,coopMwA)
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



C------------ cooperation levels ----------------------------------------------------

		subroutine calccoop(coopARwA,coopARwB,coopBRwA,coopBRwB, coopMwA)
		implicit none
      include "common_brote.h"
		double precision coopARwA,coopARwB,coopBRwA,coopBRwB,
     & 				coopMwA
		double precision p12,p21,p22, Cfun, p11A, p1A1
		integer ia, ib
		coopARwA=0.d0
		coopARwB=0.d0
		coopBRwA=0.d0
		coopBRwB=0.d0
		coopMwA=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
		  coopARwA=coopARwA
     &	            + Cfun(ia,x1AGo)*Cfun(ib,x1AGo)*iaction(2-ia,2-ib)
		  coopARwB=coopARwB
     &	            + Cfun(ia,x1AGo)*Cfun(ib,x1BGo)*iaction(2-ia,2-ib)
		  coopBRwA=coopBRwA
     &	            + Cfun(ia,x1BGo)*Cfun(ib,x1AGo)*iaction(2-ia,2-ib)
		  coopBRwB=coopBRwB
     &	            + Cfun(ia,x1BGo)*Cfun(ib,x1BGo)*iaction(2-ia,2-ib)
		  coopMwA=coopMwA
     &	            + Cfun(ia,x2oG)*Cfun(ib,x1AoG)*iaction(2-ia,2-ib)
		enddo
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






