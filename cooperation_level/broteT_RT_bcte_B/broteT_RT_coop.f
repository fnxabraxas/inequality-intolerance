c Calculates cooperation levels after invasion
c Each row ia different epsA. each column is a different epsB
c Results are neither multiplied by (1-epsA),(1-epsB) nor y,(1-y): 
c   coopAR_TOT = ( y*coopARwA ) *(1-epsA)
c   coopBR_TOT = ( y*coopBRwA + (1-y)*coopBRwB  ) *(1-epsA) 

		include "sub_common.f"

		program broteR

		implicit none
		include "common_brote.h"
		character*80 outp, outfileListEps,
     +		  outfileARtoA, outfileBRtoA, outfileBRtoB 
		integer i,j,k,stg, iepsB,iepsA, NepsAB, iepsAB
		double precision epsABini,epsABfin,
     +			coopARwA, coopBRwA, coopBRwB, epsAB, wepsAB


		print*, 'Number of strategy: '
		read(*,*) stg
		print*, 'Proportion of A (y): '
		read(*,*) yy
		print*,'epsAB: initial, final, number: '
		read(*,*) epsABini, epsABfin, NepsAB
	print*, 'Output files (4): '
	read(*,*) outfileListEps, outfileARtoA, outfileBRtoA, outfileBRtoB  
						!print*,epsAini, epsAfin, wepsA
		wepsAB=(epsABfin-epsABini)/(NepsAB-1.d0)

		call num2actionmoral(stg)
		print*,'Strategy, y:  ',stg,yy


      open(90,file=outfileListEps,status='UNKNOWN')
      epsAB=epsABini
      do iepsAB=1,NepsAB
        write(90,*) epsAB
      	epsAB=epsAB+wepsAB
      enddo
      close(90)
      
      open(91,file=outfileARtoA,status='UNKNOWN')
      open(92,file=outfileBRtoA,status='UNKNOWN')
      open(93,file=outfileBRtoB,status='UNKNOWN')
      
      epsA=epsABini
      do iepsA=1,NepsAB
        epsB=epsABini
        do iepsB=1,NepsAB
          call process(coopARwA, coopBRwA, coopBRwB)          
          write(91,'(F10.4,$)') coopARwA
          write(92,'(F10.4,$)') coopBRwA
          write(93,'(F10.4,$)') coopBRwB
          epsB=epsB+wepsAB
        enddo !iepsB
        write(91,'(A,$)')
        write(92,'(A,$)')
        write(93,'(A,$)')        
        epsA=epsA+wepsAB
      enddo !iepsA
      close(91)
      close(92)
      close(93)
      
      return
      end



		subroutine process(coopARwA, coopBRwA, coopBRwB) !(payoff11,payoff21,payoff12,payoff22)		
		implicit none
		include "common_brote.h"
		double precision coopARwA, coopBRwA, coopBRwB !payoff11,payoff21,payoff12,payoff22
		integer errorf
		double precision den,pp

					!epsA=0.3
					!yy=0.11
					!epsB=epsA

				eps=epsA
				!x1AGo=pp(yy,1.d0-yy,0)				
				call IC
				call solve_x1
				x1AoG=x(1,1,1)+x(1,2,1)
				x1AGo=x(1,1,1)+x(1,1,2)
				eps=epsA
				call calc_ppB(yy*x1AoG+(1.d0-yy)*x1BoG,1.d0,0.d0,x2oG,den)
				!if(abs(den).lt.cero) stop 'den=0 en ppal'
				if(abs(den).lt.cero) print*, 'den=0 en ppal'
				call solve_x2(errorf)
				x2Go=x(2,1,1)+x(2,1,2)
				!call calcpay(payoff11,payoff21,payoff12,payoff22)
				call calccoop(coopARwA, coopBRwA, coopBRwB)
			!print*,x1AGo,x1AoG,x2Go,x2oG
			!print*, yy, epsA, epsB	
			!print*,'Pays: ',payoff11,payoff21,payoff12,payoff22
					!stop
		return
		end


c----------- Solves ---------------------------------------------------------------



		subroutine solve_x1

      implicit none
      include "common_brote.h"
		integer a1,a2, cont, is
		double precision difx2(3),xGGold(3),xGBold(3),
     +						 		xBGold(3),xBBold(3), suma(3),
     +							 x11new,x22new, x12new,x21new, Stt,Srt,Sh

						!print*
						!print*,x(1,1,1),x(1,1,2)
						!print*,x(1,2,1),x(1,2,2)
						!print*,'---',x1AGo
						!print*,x(3,1,1),x(3,1,2)
						!print*,x(3,2,1),x(3,2,2)
						!print*,'---',x1BoG
						!print*

		cont=0
		difx2=1.d0
		x1BoG=x(3,1,1)+x(3,2,1)
		x1AoG=x(1,1,1)+x(1,2,1)
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
		  !x(1,1,2)=x1AGo-x(1,1,1)
		  !x(1,2,1)=1.d0-x1AGo-x(1,2,2)
		  !x1AoG=x(1,1,1)+x(1,2,1)

		  eps=epsB
		  x1BoG=(1.d0-yy)*Sh(x1BoG,x1BoG)+yy*Sh(x1BoG,x1AoG)

		  x(3,1,1)=0.d0
		  x(3,2,2)=1.d0-x1BoG
		  x(3,1,2)=0.d0
		  x(3,2,1)=x1BoG



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
			!if(cont.eq.50) stop

		 do is=1,3,2
		  difx2(is)=(x(is,1,1)-XGGold(is))**2.d0
     +		+(x(is,1,2)-XGBold(is))**2.d0
     +      +(x(is,2,1)-XBGold(is))**2.d0+(x(is,2,2)-XBBold(is))**2.d0
		 enddo

		enddo

		do a1=1,2
		do a2=1,2
			x(1,a1,a2)=(anint(x(1,a1,a2)*1.d0/cero))*cero
			if(x(1,a1,a2).lt.(-1.d0*cero)) then
				print*,a1,a2,x(1,a1,a2)
				print*, 'menor que 0 x1******'
				!stop
			endif
		enddo
		enddo
		x1BoG=(anint(x1BoG*1.d0/cero))*cero

		return 
		end




		subroutine solve_x2(errorf)

      implicit none
      include "common_brote.h"
		integer i,j,b1,b2,errorf
		double precision  AT(2,2),bmT(2),AR(2,2),bmR(2),A(2,2),bm(2),
     +						Ainv(2,2),xsol(2), Str

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
				print*, 'menor que 0 x2******',xsol
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
		double precision xHA, xHB, pp

		eps=epsA
		xHA=pp(1.d0,0.d0,1)
		eps=epsB
		xHB=pp(1.d0,0.d0,1)
		!print*,xH
		do i=1,3
			x(i,1,1)=xHA
			x(i,2,2)=1.d0-xHA
			x(i,1,2)=0.d0
			x(i,2,1)=0.d0
		enddo
			x(3,1,1)=xHB
			x(3,2,2)=1.d0-xHB
		!do i=1,3
		!	x(i,1,1)=x1AGo
		!	x(i,2,2)=1.d0-x1AGo
		!	x(i,1,2)=0.d0
		!	x(i,2,1)=0.d0
		!enddo
		return
		end



C------------ cooperation levels ----------------------------------------------------

		subroutine calccoop(coopARwA, coopBRwA, coopBRwB)
		implicit none
      include "common_brote.h"
		double precision  coopARwA, coopBRwA, coopBRwB
		double precision p12,p21,p22, Cfun, p11A, p1A1
		integer ia, ib
		coopARwA=0.d0
		coopBRwA=0.d0
		coopBRwB=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
		  coopARwA=coopARwA
     &	            + Cfun(ia,x1AGo)*Cfun(ib,x1AGo)*iaction(2-ia,2-ib) 
		  coopBRwA=coopBRwA
     &	            + Cfun(ia,x1BoG)*Cfun(ib,x1AoG)*iaction(2-ia,2-ib)
		  coopBRwB=coopBRwB
     &	            + Cfun(ia,x1BoG)*Cfun(ib,x1BoG)*iaction(2-ia,2-ib)
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
		double precision p11,p12,p21,p22, Cfun, p11A, p1A1
		integer ia, ib

			
		p11A=0.d0
		p1A1=0.d0
		p12=0.d0
		p21=0.d0
		p22=0.d0
		do ia=1,0,-1
		do ib=1,0,-1
			p1A1=p1A1+ ( Cfun(ia,x1AGo)*yy*Cfun(ib,x1AGo)*(1.d0-epsA)  
     &	  							)*iaction(2-ia,2-ib)        ! 1A le da a 1
			p11A=p11A+ ( Cfun(ia,x1AGo)*yy*Cfun(ib,x1AGo)*(1.d0-epsA)
     &	 + (1.d0-yy)*Cfun(ia,x1BoG)*Cfun(ib,x1AoG)*(1.d0-epsB)	
     &								)*iaction(2-ia,2-ib)  	! 1 le da a 1A
			p22=p22+ Cfun(ia,x2oG)*Cfun(ib,x2oG)*iaction(2-ia,2-ib) 	! 2 le da a 2
			p21=p21+ ( Cfun(ia,x2oG)*Cfun(ib,x1AoG)*yy
     &   + Cfun(ia,x2oG)*Cfun(ib,x1BoG)*(1.d0-yy)  	
     &								)*iaction(2-ia,2-ib)	! 2 le da a 1
			p12=p12+ ( Cfun(ia,x1AGo)*yy*Cfun(ib,x2Go)*(1.d0-epsA)   
     &	 + (1.d0-yy)*Cfun(ia,x1BoG)*Cfun(ib,x2oG)*(1.d0-epsB)  
     &								)*iaction(2-ia,2-ib) 	! 1 le da a 2
		enddo
		enddo
		p21=p21 *(1.d0-epsA) 
		p22=p22 *(1.d0-epsA)

				!print*,p22,p21,p12,p11A,p11B
			!print*,'p: ', p11A,p1A1,x1AGo,x1AoG

			payoff11=b*p11A -c*p1A1   ! de 1 (con 1)
			payoff21=(b*p12-c*p21)   ! de 2 (con 1)
			payoff12=(b*p21-c*p12)   ! de 1 (con 2)
			payoff22=((b-c)*p22  ) 	 ! de 2 (con 2)


		return
		end



