

		integer, parameter :: Nbmax=1000, Nepsmax=1000, Nymax=1000
		double precision, parameter :: c=1.d0, cero=1.d-15,  ! d-15 normal
     +					dt=1.d-1, xtol=1.d-10, epstol=1.d-4  ! xtol: d-12 normal

		integer  iaction(2,2),imoral(2,2,2), Nb, Neps, NepsB,Ny
		double precision eps, x1AGo, x1AoG,x1BGo,x1BoG, x2oG, x2Go,b,
     +			 yy, bvect(11), yvect(11),
     +				x(3,2,2), epsB, epsA, epsBvect(Nbmax)

	   common /ppal/ eps,epsB,epsA,yy,bvect,b,
     +					yvect,epsBvect,Nb,Neps,Ny
	   common /xx/x1AGo,x1AoG,x1BGo,x1BoG,x2oG,x2Go,iaction,imoral,x
