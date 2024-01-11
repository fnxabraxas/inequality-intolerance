

		integer, parameter :: Nbmax=1000, Nepsmax=1000, Nymax=1000
		double precision, parameter :: c=1.d0, cero=1.d-12,  ! d-15 normal
     +					dt=1.d-1, xtol=1.d-12, epstol=1.d-4

		integer  iaction(2,2),imoral(2,2,2), Nb, Neps, NepsA,Ny
		double precision eps, x1AGo, x1AoG,x1BGo,x1BoG, x2oG, x2Go,b,
     +			 yy, bvect(10), yvect(10),
     +				x(3,2,2), epsB, epsA, epsAvect(Nbmax)

	   common /ppal/ eps,epsB,epsA,yy,bvect,b,
     +					yvect,epsAvect,Nb,Neps,Ny
	   common /xx/x1AGo,x1AoG,x1BGo,x1BoG,x2oG,x2Go,iaction,imoral,x
