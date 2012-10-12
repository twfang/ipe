!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
!      module sub_efield
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use physconst,     only: pi
c     use abortutils,    only: endrun
c     use cam_logfile,   only: iulog
   
!      implicit none

!nm20121003:module parameters are separated into efield.f!
!      end module sub_efield
!
!      
! Purpose: 
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
! Method:
!
! To use, first call subroutine ReadCoef once.
! Next, call SetModel with the specified input parameters.
! The function EpotVal(gLAT,gMLT) can then be used repeatively to get the
! electric potential at the desired location in geomagnetic coordinates.
! Subroutines to calculate the electric potentials from the Weimer '96 model of
! the polar cap ionospheric electric potentials.
!
!
! Author: A. Maute Dec 2003  
! This code is protected by copyright and is
! distributed for research or educational use only.
! Commerical use without written permission from Dan Weimer/MRC is prohibited.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!==================================================================

	FUNCTION EpotVal(gLAT,gMLT)
!
!-----------------------------------------------------------------------
! Return the value of the electric potential in kV at
! corrected geomagnetic coordinates gLAT (degrees) and gMLT (hours).
!
! Must first call ReadCoef and SetModel to set up the model coeficients for
! the desired values of Bt, IMF clock angle, Dipole tilt angle, and SW Vel.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real EpotVal
!
!-------------------------------Commons---------------------------------
!
	INTEGER ML,MM
	REAL Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL gLAT,gMLT
!
!---------------------------Local variables-----------------------------
!
        integer limit,l,m

	Real Theta,Phi,Z,ct,Phim
        real r
	REAL Plm(0:20,0:20)
!
!-----------------------------------------------------------------------
!
	r=90.-gLAT
	IF(r .LT. 45.)THEN
	  Theta=r*pi/45.
          Phi=gMLT*pi/12.
	  Z=Coef(0,0,0)
	  ct=COS(Theta)
	  CALL Legendre(ct,ML,MM,Plm)
	  DO l=1,ML
	    Z=Z + Coef(0,l,0)*Plm(l,0)
	    IF(l.LT.MM)THEN
	      limit=l
	    ELSE
	      limit=MM
	    ENDIF
	    DO m=1,limit
	      phim=phi*m
              Z=Z + Coef(0,l,m)*Plm(l,m)*COS(phim) +  
     &	   Coef(1,l,m)*Plm(l,m)*SIN(phim)
	    ENDDO
	  ENDDO
	ELSE
	  Z=0.
	ENDIF
	EpotVal=Z
	RETURN
	END FUNCTION EpotVal

!================================================================================================

	SUBROUTINE ReadCoef (wei96_file)
!
!-----------------------------------------------------------------------
!
! Read in the data file with the model coefficients
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!
! NCAR addition (Jan 97):  initialize constants used in GECMP
!-----------------------------------------------------------------------
!
c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use ioFileMod,     only : getfil
c     use units,         only : getunit, freeunit
c     use abortutils,    only : endrun
c     use cam_logfile,   only : iulog
      implicit none 
!
!-------------------------------Commons---------------------------------
!
      real alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!            ALAMN = Absolute min latitude (deg) of model
!            ALAMX = Absolute max latitude (deg) for normal gradient calc.
!            STPD  = Angular dist (deg) of step @ 300km above earth (r=6371km)
!            STP2  = Denominator in gradient calc

!
!------------------------------Arguments--------------------------------
!
      character(len=*), intent(in) :: wei96_file
!
!-----------------------------Parameters------------------------------
!
      real d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847 ,  
     &            R2D = 57.2957795130823208767981548147)
!
!---------------------------Local variables-----------------------------
!
      INTEGER udat,unit,ios
      integer ll,mm,k,m,klimit,kk,nn,ii,i,n,ilimit,mlimit,l

      REAL C(0:3)
      real stpr, step

      CHARACTER*15 skip

      INTEGER MaxL,MaxM,MaxN,iulog
      REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
      COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

      character(len=256) :: locfn
!
!-----------------------------------------------------------------------
      iulog=14
      STEP = 10.
      STPR = STEP/6671.
      STPD = STPR*R2D
      STP2 = 2.*STEP
      CSTP = COS (STPR)
      SSTP = SQRT (1. - CSTP*CSTP)
      ALAMN = 45.
      ALAMX = 90. - STPD
      ALAMR = ALAMN*D2R
!          End NCAR addition
! 
!  get coeff_file  
c     unit= getunit()
      unit= 200
c     write(iulog,*) 'Weimer: getting file ',trim(wei96_file),
c    &' unit ',unit
c     call getfil( wei96_file, locfn, 0 )
      locfn= wei96_file
!      
c     write(iulog,*) 'Weimer: opening file ',trim(locfn),
c    &' unit ',unit	
      OPEN(unit=unit,file=trim(locfn),  
     &     status = 'old',iostat = ios)
c     if(ios.gt.0) then
c      write(iulog,*) 'Weimer: error in opening wei96.cofcnts',
c    &' unit ',unit
c       call endrun
c     endif
      
  900 FORMAT(A15)
 1000 FORMAT(3I8)
 2000 FORMAT(3I2)
 3000 FORMAT(2I2,4E15.6)

!     READ(udat,900) skip
c     write(iulog,*) 'Weimer: reading file ',trim(locfn),
c    &' unit ',unit	
      READ(unit,1000,iostat = ios) MaxL,MaxM,MaxN
c     if(ios.gt.0) then
c     write(iulog,*) 
c    &'ReadCoef: error in reading wei96.cofcnts file',
c    &' unit ',unit	
c       call endrun
c     endif
      DO l=0,MaxL
        IF(l.LT.MaxM)THEN
          mlimit=l
        ELSE
          mlimit=MaxM
        ENDIF
        DO m=0,mlimit
          IF(m.LT.1)THEN
            klimit=0
          ELSE
            klimit=1
          ENDIF
          DO k=0,klimit
            READ(unit,2000,iostat = ios) ll,mm,kk
c           if(ios.gt.0) then
c     	      write(iulog,*) 
c    &'ReadCoef: error in reading wei96.cofcnts file',' unit ',
c    &unit	
c             call endrun
c           endif
c           IF(ll.NE.l .OR. mm.NE.m .OR. kk.NE.k)THEN
c             WRITE(IULOG,*)'Data File Format Error'
c             CALL ENDRUN
c           ENDIF
            DO n=0,MaxN
              IF(n.LT.1)THEN
        	ilimit=0
              ELSE
        	ilimit=1
              ENDIF
              DO i=0,ilimit
        	READ(unit,3000,iostat = ios) nn,ii,C
c               if(ios.gt.0) then
c     	          write(iulog,*) 'ReadCoef: error in reading',  
c    &                 ' wei96.cofcnts file',' unit ',unit	
c                 call endrun
c               endif
c       	IF(nn.NE.n .OR. ii.NE.i)THEN
c       	  WRITE(IULOG,*)'Data File Format Error'
c       	  CALL ENDRUN
c       	ENDIF
        	Cn(0,i,n,k,l,m)=C(0)
        	Cn(1,i,n,k,l,m)=C(1)
        	Cn(2,i,n,k,l,m)=C(2)
        	Cn(3,i,n,k,l,m)=C(3)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!      
      close(unit)
c     call freeunit(unit)
!    
      RETURN
      END SUBROUTINE ReadCoef

!================================================================================================

	FUNCTION FSVal(omega,MaxN,FSC)
!
!-----------------------------------------------------------------------
! Evaluate a  Sine/Cosine Fourier series for N terms up to MaxN
! at angle omega, given the coefficients in FSC
!
!*********************** Copyright 1996, Dan Weimer/MRC ***************
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real FSVal
!
!------------------------------Arguments--------------------------------
!
	INTEGER MaxN
	REAL omega,FSC(0:1,0:*)
!
!---------------------------Local variables-----------------------------
!
	INTEGER n
	REAL Y,theta
!
!-----------------------------------------------------------------------
!
	Y=0.
	DO n=0,MaxN
	  theta=omega*n
	  Y=Y + FSC(0,n)*COS(theta) + FSC(1,n)*SIN(theta)
	ENDDO
	FSVal=Y
	RETURN
	END FUNCTION FSVal

!================================================================================================

	SUBROUTINE SetModel(angle,Bt,Tilt,SWVel)
!
!-----------------------------------------------------------------------
! Calculate the complete set of spherical harmonic coefficients,
! given an arbitrary IMF angle (degrees from northward toward +Y),
! magnitude Bt (nT), dipole tilt angle (degrees),
! and solar wind velocity (km/sec).
! Returns the Coef in the common block SetCoef.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-------------------------------Commons---------------------------------
!
	INTEGER MaxL,MaxM,MaxN
	REAL Cn( 0:3 , 0:1 , 0:4 , 0:1 , 0:8 , 0:3 )
	COMMON /AllCoefs/Cn,MaxL,MaxM,MaxN

	INTEGER ML,MM
	REAL Coef(0:1,0:8,0:3),pi
	COMMON/SetCoef/Coef,pi,ML,MM
!
!------------------------------Arguments--------------------------------
!
	REAL angle,Bt,Tilt,SWVel
!
!---------------------------Local variables-----------------------------
!
        integer n, k, ilimit, i, klimit, l, m, mlimit
	REAL FSC(0:1,0:4), fsval, omega, sintilt
!
!-----------------------------------------------------------------------
!
	pi=2.*ASIN(1.)
	ML=MaxL
	MM=MaxM
	SinTilt=SIN(Tilt*pi/180.)
!	SinTilt=SIND(Tilt)

	omega=angle*pi/180.

        fsc(1,0) = 0.
	DO l=0,MaxL
	  IF(l.LT.MaxM)THEN
	    mlimit=l
	  ELSE
	    mlimit=MaxM
	  ENDIF
	  DO m=0,mlimit
	    IF(m.LT.1)THEN
	      klimit=0
	    ELSE
	      klimit=1
	    ENDIF
	    DO k=0,klimit
! Retrieve the regression coefficients and evaluate the function
! as a function of Bt,Tilt,and SWVel to get each Fourier coefficient.
	      DO n=0,MaxN
	        IF(n.LT.1)THEN
	          ilimit=0
	        ELSE
	          ilimit=1
	        ENDIF
		DO i=0,ilimit
		  FSC(i,n)=Cn(0,i,n,k,l,m) + Bt*Cn(1,i,n,k,l,m) +  
     &	   SinTilt*Cn(2,i,n,k,l,m) + SWVel*Cn(3,i,n,k,l,m)
		ENDDO
	      ENDDO
! Next evaluate the Fourier series as a function of angle.
      	      Coef(k,l,m)=FSVal(omega,MaxN,FSC)
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE SetModel

!================================================================================================

	SUBROUTINE LEGENDRE(x,lmax,mmax,Plm)
!
!-----------------------------------------------------------------------
! compute Associate Legendre Function P_l^m(x)
! for all l up to lmax and all m up to mmax.
! returns results in array Plm
! if X is out of range ( abs(x)>1 ) then value is returned as if x=1.
!
!*********************** Copyright 1996, Dan Weimer/MRC ***********************
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
c       use cam_logfile,  only : iulog

        implicit none 
!
!------------------------------Arguments--------------------------------
!
        integer lmax, mmax
	real x, Plm(0:20,0:20)
!
!---------------------------Local variables-----------------------------
!
        integer m, lm2, l, iulog
        real xx, fact
        iulog=14
!
!-----------------------------------------------------------------------
!
	  DO l=0,20
	    DO m=0,20
		Plm(l,m)=0.
	    ENDDO
	  ENDDO
	xx=MIN(x,1.)
	xx=MAX(xx,-1.)
c         IF(lmax .LT. 0 .OR. mmax .LT. 0 .OR. mmax .GT. lmax )THEN
c        write(iulog,*)'Bad arguments to Legendre'
c        RETURN
c        ENDIF
! First calculate all Pl0 for l=0 to l
	Plm(0,0)=1.
	IF(lmax.GT.0)Plm(1,0)=xx
	IF (lmax .GT. 1 )THEN
	  DO L=2,lmax
	    Plm(L,0)=( (2.*L-1)*xx*Plm(L-1,0) - 
     &(L-1)*Plm(L-2,0) )/L
	  ENDDO
	ENDIF
	IF (mmax .EQ. 0 )RETURN
	fact=SQRT( (1.-xx)*(1.+xx) )
	DO M=1,mmax
	  DO L=m,lmax
	    lm2=MAX(L-2,0)
	    Plm(L,M)=Plm(lm2,M) - ( 2*L-1)*fact*Plm(L-1,M-1)
	  ENDDO
	ENDDO
	RETURN
	END SUBROUTINE LEGENDRE

!================================================================================================

!*********************** Copyright 1996, Dan Weimer/MRC ***********************

!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following routines (translib.for) were added to return the dipole tilt. C
!  GET_TILT was initially a procedure (TRANS), here it has been changed into   C
!  a function which returns the dipole tilt.                                   C
! Barbara Emery (emery@ncar.ucar.edu) and William Golesorkhi, HAO/NCAR (3/96)  C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! COORDINATE TRANSFORMATION UTILITIES
!**********************************************************************        
	FUNCTION GET_TILT(YEAR,MONTH,DAY,HOUR)
!
!-----------------------------------------------------------------------
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The following line initially was:                                           C
!       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)                            C
!  It has been changed to return the dipole tilt from this function call.      C
!CC NCAR MODIFIED (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         
!      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!	   K=10    "      GEO to SM 
!	   K=11    "      MAG to SM 
!
!      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
!      FILE UNIT=IDBUG
!
!      The formal names of the coordinate systems are:
!	GSE - Geocentric Solar Ecliptic
!	GEO - Geographic
!	MAG - Geomagnetic
!	GSM - Geocentric Solar Magnetospheric
!	SM  - Solar Magnetic
!	
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!
!      Program author:  D. R. Weimer
!
!      Some of this code has been copied from subroutines which had been
!      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
!      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
!
!      The formulas for the calculation of Greenwich mean sidereal time (GMST)
!      and the sun's location are from "Almanac for Computers 1990",
!      U.S. Naval Observatory.
!
!-----------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value--------------------------
!
        real get_tilt
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

        real epoch, th0, ph0, dipole
        COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
c       DATA EPOCH,TH0,PH0,DIPOLE/1980.,11.19,-70.76,.30574/
!
!------------------------------Arguments--------------------------------
!
	INTEGER YEAR, MONTH, DAY
	REAL HOUR
!
!-----------------------------Parameters------------------------------
!
	INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
	INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
	INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
	INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
	INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

	PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
	PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
	PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
	PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
	PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
	PARAMETER (MAGSM =11,SMMAG =-11)
!
!---------------------------Local variables-----------------------------
!
        integer IDBUG
        integer j, k, jd, iyr, i, mjd

        REAL UT, T0, GMSTD, GMSTH, ECLIP, MA, LAMD, SUNLON, pi
        real b32, b33, b3
!
!-------------------------External Functions----------------------------
!
        integer julday
        external julday
!
!-----------------------------------------------------------------------
!
        EPOCH=1980.
        TH0=11.19
        PH0=-70.76
        DIPOLE=.30574
	pi=2.*ASIN(1.)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IDBUG=0 to prevent printing data to the screen or writing data to a file.    C
        IDBUG = 0
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	IF(YEAR.LT.1900)THEN
	  IYR=1900+YEAR
	ELSE
	  IYR=YEAR
	ENDIF
	UT=HOUR
	JD=JULDAY(MONTH,DAY,IYR)
	MJD=JD-2400001
c       T0=(real(MJD,r8)-51544.5)/36525.0
	T0=(float(MJD)-51544.5)/36525.0
	GMSTD=100.4606184 +36000.770*T0 +3.87933E-4*T0*T0 + 
     &       15.0410686*UT
	CALL ADJUST(GMSTD)
	GMSTH=GMSTD*24./360.
	ECLIP=23.439 - 0.013*T0
        MA=357.528 + 35999.050*T0 + 0.041066678*UT
        CALL ADJUST(MA)
        LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
        CALL ADJUST(LAMD)
        SUNLON=LAMD + (1.915-0.0048*T0)*SIN(MA*pi/180.) + 0.020* 
     &     SIN(2.*MA*pi/180.)
        CALL ADJUST(SUNLON)
c         IF(IDBUG.NE.0)THEN
c         WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
c         WRITE(IDBUG,*) 'MJD=',MJD
c         WRITE(IDBUG,*) 'T0=',T0
c         WRITE(IDBUG,*) 'GMSTH=',GMSTH
c         WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
c         WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
c         WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
c         WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
c         ENDIF

	CX(1)= GMSTD
	CX(2) = ECLIP
	CX(3) = SUNLON
	CX(4) = TH0
	CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

	DO I=1,5
	  ST(I) = SIN(CX(I)*pi/180.)
	  CT(I) = COS(CX(I)*pi/180.)
	ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.         
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
!         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0.         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0.         
      AM(3,1,GEIGEO) = 0.         
      AM(3,2,GEIGEO) = 0.         
      AM(3,3,GEIGEO) = 1.         
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) + 
     &AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
!         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) + 
     &AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0.) B3 = -B3  
!         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) + 
     &AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + 
     &AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) + 
     &AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + 
     &AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO 
!
	ST(6) = AM(3,1,GSEMAG)        
	CT(6) = SQRT(1.-ST(6)*ST(6))      
	CX(6) = ASIN(ST(6)*pi/180.)  

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0.
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0.
        AM(2,2,GSMSM) = 1.
        AM(2,3,GSMSM) = 0.
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0.
        AM(3,3,GSMSM) = CT(6)  
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) +  
     &AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  
     &AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) +  
     &  AM(I,2,GSMSM)*AM(J,2,GSMMAG) + 
     &AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
      
!
      CX(7)=ATAN2( AM(2,1,11) , AM(1,1,11) )
      
      CX(7)=CX(7)*180./pi
      CX(8)=ASIN( AM(3,1,1)*pi/180. )
      CX(9)=ATAN2( AM(2,1,1) , AM(1,1,1) )
      CX(9)=CX(9)*180./pi

      IF(IDBUG.NE.0)THEN
c     WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
c     WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',
c    &CX(7)
c     WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
c     WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
c        WRITE(IDBUG,1001) K
         DO I=1,3
c          WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  The next line was added to return the dipole tilt from this function call.  C
      GET_TILT = CX(6)
!CC NCAR MODIFICATION (3/96) CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END FUNCTION GET_TILT

!======================================================================

      SUBROUTINE ROTATE (X,Y,Z,I)       
!
!-----------------------------------------------------------------------
!     THIS SUBROUTINE APPLIES TO THE VECTOR (X,Y,Z) THE ITH ROTATION  
!     MATRIX AM(N,M,I) GENERATED BY SUBROUTINE TRANS
!     IF I IS NEGATIVE, THEN THE INVERSE ROTATION IS APPLIED
!-----------------------------------------------------------------------
!
c     use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!------------------------------Arguments--------------------------------
!
      integer i
      REAL X,Y,Z
!
!---------------------------Local variables-----------------------------
!
      REAL A(3)
!
!-----------------------------------------------------------------------
!
      A(1)=X
      A(2)=Y
      A(3)=Z
      CALL ROTATEV(A,A,I)
      X=A(1)
      Y=A(2)
      Z=A(3)
    
      RETURN        
      END SUBROUTINE ROTATE

!======================================================================

      SUBROUTINE ROTATEV (A,B,I)       
!         
!-----------------------------------------------------------------------
!     THIS SUBROUTINE APPLIES TO THE VECTOR A(3) THE ITH ROTATION  
!     MATRIX AM(N,M,I) GENERATED BY SUBROUTINE TRANS
!     AND OUTPUTS THE CONVERTED VECTOR B(3), WITH NO CHANGE TO A.
!     IF I IS NEGATIVE, THEN THE INVERSE ROTATION IS APPLIED
!-----------------------------------------------------------------------
!
c     use shr_kind_mod, only: r8 => shr_kind_r8
c     use cam_logfile,  only : iulog
c     use abortutils,   only : endrun

      implicit none 
!
!-------------------------------Commons---------------------------------
!
      real cx, st, ct, am
      COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
      integer i
      REAL A(3),B(3)
!
!---------------------------Local variables-----------------------------
!
      integer id, j, iulog
      real xa, ya, za
      iulog=14
!
!-----------------------------------------------------------------------
!
c     IF(I.EQ.0 .OR. IABS(I).GT.11)THEN
c     WRITE(IULOG,*)'ROTATEV CALLED WITH UNDEFINED TRANSFORMATION'
c     CALL ENDRUN
c     ENDIF

      XA = A(1)
      YA = A(2)
      ZA = A(3)
      IF(I.GT.0)THEN
	ID=I
        DO J=1,3
          B(J) = XA*AM(J,1,ID) + YA*AM(J,2,ID) + ZA*AM(J,3,ID)
        ENDDO
      ELSE
	ID=-I
        DO J=1,3
          B(J) = XA*AM(1,J,ID) + YA*AM(2,J,ID) + ZA*AM(3,J,ID)
        ENDDO
      ENDIF
      RETURN        
      END SUBROUTINE ROTATEV

!================================================================================================

	SUBROUTINE FROMCART(R,LAT,LONG,POS)
!
!-----------------------------------------------------------------------
! CONVERT CARTESIAN COORDINATES POS(3)
! TO SPHERICAL COORDINATES R, LATITUDE, AND LONGITUDE (DEGREES)
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
	REAL R, LAT, LONG, POS(3)
!
!---------------------------Local variables-----------------------------
!
        real pi
!
!-----------------------------------------------------------------------
!
	pi=2.*ASIN(1.)
	R=SQRT(POS(1)*POS(1) + POS(2)*POS(2) + POS(3)*POS(3))
	IF(R.EQ.0.)THEN
	  LAT=0.
	  LONG=0.
	ELSE
	  LAT=ASIN(POS(3)*pi/180./R)
	  LONG=ATAN2(POS(2),POS(1))
	  LONG=LONG*180./pi
	ENDIF
	RETURN
	END SUBROUTINE FROMCART

!================================================================================================

	SUBROUTINE TOCART(R,LAT,LONG,POS)
!
!-----------------------------------------------------------------------
! CONVERT SPHERICAL COORDINATES R, LATITUDE, AND LONGITUDE (DEGREES)
! TO CARTESIAN COORDINATES POS(3)
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
	REAL R, LAT, LONG, POS(3)
!
!---------------------------Local variables-----------------------------
!
        real pi, stc, ctc, sf, cf
!
!-----------------------------------------------------------------------
!
	pi=2.*ASIN(1.)
        STC = SIN(LAT*pi/180.)    
        CTC = COS(LAT*pi/180.)    
        SF = SIN(LONG*pi/180.)     
        CF = COS(LONG*pi/180.)     
        POS(1) = R*CTC*CF        
        POS(2) = R*CTC*SF        
        POS(3) = R*STC
	RETURN
	END SUBROUTINE TOCART

!================================================================================================

	SUBROUTINE ADJUST(ANGLE)
!
!-----------------------------------------------------------------------
!	ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!------------------------------Arguments--------------------------------
!
        real angle
!
!-----------------------------------------------------------------------
!
 10	CONTINUE
	IF(ANGLE.LT.0.)THEN
	  ANGLE=ANGLE+360.
	  GOTO 10
	ENDIF
 20	CONTINUE
 	IF(ANGLE.GE.360.)THEN
	  ANGLE=ANGLE-360.
	  GOTO 20
	ENDIF
	RETURN
	END SUBROUTINE ADJUST

!================================================================================================

      INTEGER FUNCTION JULDAY(MM,ID,IYYY)
!
!-----------------------------------------------------------------------
!
c     use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!------------------------------Arguments--------------------------------
!
      integer mm, id, iyyy
!
!-----------------------------Parameters------------------------------
!
      integer igreg
      PARAMETER (IGREG=15+31*(10+12*1582))
!
!---------------------------Local variables-----------------------------
!
      integer ja, jm, jy
!
!-----------------------------------------------------------------------
!
      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
      END FUNCTION JULDAY

!================================================================================================

	FUNCTION MLT(MagLong)
!
!-----------------------------------------------------------------------
! given magnetic longitude in degrees, return Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real mlt
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
	COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)

!
!------------------------------Arguments--------------------------------
!
	REAL MagLong
!
!---------------------------Local variables-----------------------------
!
	REAL angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
!       MLT=ADJUST(Maglong+RotAngle+180.)/15.
        angle = Maglong+RotAngle+180.
        call adjust(angle)
        mlt = angle/15.
	RETURN
	END FUNCTION MLT

!================================================================================================

	FUNCTION MagLong(MLT)
!
!-----------------------------------------------------------------------
! return magnetic longitude in degrees, given Magnetic Local Time
! assuming that TRANS has been called with the date & time to calculate
! the rotation matrices.
!
! btf 11/06/03:
! Call sub adjust instead of referencing it as a function
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-----------------------------Return Value------------------------------
!
        real MagLong
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	REAL MLT
!
!---------------------------Local variables-----------------------------
!
	REAL angle, rotangle
!
!-----------------------------------------------------------------------
!
	RotAngle=CX(7)
	angle=MLT*15.-RotAngle-180.
!       MagLong=ADJUST(angle)
        call adjust(angle)
        MagLong = angle
	RETURN
	END FUNCTION MagLong

!================================================================================================

	SUBROUTINE SunLoc(SunLat,SunLong)
!
!-----------------------------------------------------------------------
! Return latitude and longitude of sub-solar point.
! Assumes that TRANS has previously been called with the
! date & time to calculate the rotation matrices.
!-----------------------------------------------------------------------
!
c       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none 
!
!-------------------------------Commons---------------------------------
!
        real cx, st, ct, am
        COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!
!------------------------------Arguments--------------------------------
!
	Real SunLat,SunLong
!
!-----------------------------------------------------------------------
!
	SunLong=CX(9)
	SunLat=CX(8)
	RETURN
	END SUBROUTINE SunLoc

!================================================================================================

      SUBROUTINE GECMP (AMLA,RMLT,ET,EP)
!
!-----------------------------------------------------------------------
!          Get Electric field components for the Weimer electrostatic
!          potential model.  Before use, first load coefficients (CALL
!          READCOEF) and initialize model conditions (CALL SETMODEL).
!
!          INPUTS:
!            AMLA = Absolute value of magnetic latitude (deg)
!            RMLT = Magnetic local time (hours).
!          RETURNS:
!            ET = Etheta (magnetic equatorward*) E field component (V/m)
!            EP = Ephi   (magnetic eastward)     E field component (V/m)
!
!          * ET direction is along the magnetic meridian away from the
!            current hemisphere; i.e., when ET > 0, the direction is
!              southward when RMLA > 0
!              northward when RMLA < 0
!
!          NCAR addition (Jan 97).  R.Barnes
!-----------------------------------------------------------------------
!
c     use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!
!-------------------------------Commons---------------------------------
!
!          CECMP contains constants initialized in READCOEF
      real alamn, alamx, alamr, stpd, stp2, cstp, sstp
      COMMON /CECMP/ ALAMN,ALAMX,ALAMR,STPD,STP2,CSTP,SSTP
!
!------------------------------Arguments--------------------------------
!
      real amla, rmlt, et, ep
!
!-----------------------------Parameters------------------------------
!
      real d2r, r2d
      PARAMETER ( D2R =  0.0174532925199432957692369076847 , 
     &           R2D = 57.2957795130823208767981548147)
!
!---------------------------Local variables-----------------------------
!
      real p1, p2
      real xmlt, xmlt1, kpol, dphi, amla1
!
!-------------------------External Functions----------------------------
!
      real epotval
      external epotval
!
!-----------------------------------------------------------------------
!
      ET = -99999.
      EP = -99999.
      IF (AMLA .LT. 0.) GO TO 100

!          Calculate -(latitude gradient) by stepping 10 km along the
!          meridian in each direction (flipping coordinates when going
!          over pole to keep lat <= 90).
      KPOL  = 0
      XMLT  = RMLT
   10 XMLT1 = XMLT
      AMLA1 = AMLA + STPD
      IF (AMLA1 .GT. 90.) THEN
	AMLA1 = 180. - AMLA1
	XMLT1 = XMLT1 + 12.
      ENDIF
      P1 = EPOTVAL (AMLA1    ,XMLT1)
      P2 = EPOTVAL (AMLA-STPD,XMLT )
      IF (KPOL .EQ. 1) GO TO 20
      ET = (P1 - P2) / STP2

!          Calculate -(lon gradient).  For most latitudes, step along a
!          great circle.  However, limit minimum latitude to the model
!          minimum (distorting the path onto a latitude line).  Also,
!          avoid a divide by zero at the pole avoid by using Art's trick
!          where Ephi(90,lon) = Etheta(90,lon+90)
      IF (AMLA .LT. ALAMX) THEN
	AMLA1 = MAX (ASIN(SIN(AMLA*D2R)*CSTP) , ALAMR)
	DPHI  = ASIN (SSTP/SIN(AMLA1))*R2D
	AMLA1 = AMLA1*R2D
	P1 = EPOTVAL (AMLA1,XMLT+DPHI)
	P2 = EPOTVAL (AMLA1,XMLT-DPHI)
      ELSE
	AMLA = 90.
	XMLT = XMLT + 6.
	KPOL = 1
	GO TO 10
      ENDIF
   20 EP = (P2 - P1) / STP2
      IF (KPOL .EQ. 1) EP = -EP

!          Below model minimum lat, the potential is value at min lat
      IF (AMLA .LT. ALAMN) THEN
	ET = 0.
	EP = EP * COS(ALAMR)/COS(AMLA*D2R)
      ENDIF

  100 RETURN
      END SUBROUTINE GECMP

!=====================================================================
      subroutine svdcmp( a, m, n, mp, np, w, v )
!------------------------------------------------------------------------- 
! purpose: singular value decomposition
!
! method:
! given a matrix a(1:m,1:n), with physical dimensions mp by np,
! this routine computes its singular value decomposition,
! the matrix u replaces a on output. the
! diagonal matrix of singular values w is output as a vector
! w(1:n). the matrix v (not the transpose v^t) is output as
! v(1:n,1:n).
!
! author: a. maute dec 2003      
! (* copyright (c) 1985 numerical recipes software -- svdcmp *!
! from numerical recipes 1986 pp. 60 or can be find on web-sites
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      integer, intent(in)     :: m
      integer, intent(in)     :: n
      integer, intent(in)     :: mp
      integer, intent(in)     :: np
      real, intent(inout) :: a(mp,np)
      real, intent(out)   :: v(np,np)
      real, intent(out)   :: w(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, its, j, k, l, nm
      real :: anorm
      real  :: c
      real  :: f
      real  :: g
      real  :: h
      real  :: s
      real  :: scale
      real  :: x, y, z
      real  :: rv1(nmax)
      logical  :: cnd1
      logical  :: cnd2

      g     = 0.0
      scale = 0.0
      anorm = 0.0

      do i = 1,n  !loop1
        l = i + 1
        rv1(i) = scale*g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m ) then
          do k = i,m
            scale = scale + abs(a(k,i))
          end do
          if( scale /= 0.0 ) then
            do k = i,m
              a(k,i) = a(k,i)/scale
              s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            if( i /= n ) then
              do j = l,n
                s = 0.0
                do k = i,m
                  s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k = i,m
                  a(k,j) = a(k,j) + f*a(k,i)
                end do
              end do
            end if
            do k = i,m
              a(k,i) = scale*a(k,i)
            end do
          endif
        endif
        w(i) = scale *g
        g     = 0.0
        s     = 0.0
        scale = 0.0
        if( i <= m .and. i /= n ) then
          do k = l,n
            scale = scale + abs(a(i,k))
          end do
          if( scale /= 0.0 ) then
            do k = l,n
              a(i,k) = a(i,k)/scale
              s      = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l,n
              rv1(k) = a(i,k)/h
            end do
            if( i /= m ) then
              do j = l,m
                s = 0.0
                do k = l,n
                  s = s + a(j,k)*a(i,k)
                end do
                do k = l,n
                  a(j,k) = a(j,k) + s*rv1(k)
                end do
              end do
            end if
            do k = l,n
              a(i,k) = scale*a(i,k)
            end do
          end if
        end if
        anorm = max( anorm,(abs(w(i)) + abs(rv1(i))) )
      enddo !loop1

      do i = n,1,-1
        if( i < n ) then
          if( g /= 0.0 ) then
            do j = l,n
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l,n
              s = 0.0
              do k = l,n
                s = s + a(i,k)*v(k,j)
              end do
              do k = l,n
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do
          end if
          do j = l,n
            v(i,j) = 0.0
            v(j,i) = 0.0
          end do
        end if
        v(i,i) = 1.0
        g = rv1(i)
        l = i
      end do

      do i = n,1,-1
        l = i + 1
        g = w(i)
        if( i < n ) then
          do j = l,n
            a(i,j) = 0.0
          end do
        end if
        if( g /= 0.0  ) then
          g = 1./g
          if( i /= n ) then
            do j = l,n
              s = 0.0
              do k = l,m
                s = s + a(k,i)*a(k,j)
              end do
              f = (s/a(i,i))*g
              do k = i,m
                a(k,j) = a(k,j) + f*a(k,i)
              end do
            end do
          end if
          do j = i,m
            a(j,i) = a(j,i)*g
          end do
        else
          do j = i,m
            a(j,i) = 0.0
          end do
        end if
        a(i,i) = a(i,i) + 1.0
      end do

      do k = n,1,-1
        do its = 1,30 !loop2
          do l = k,1,-1
            nm = l - 1
            cnd1 = abs( rv1(l) ) + anorm == anorm
            if( cnd1 ) then
              cnd2 = .false.
              exit
            end if
            cnd2 = abs( w(nm) ) + anorm == anorm
            if( cnd2 ) then
              cnd1 = .true.
              exit
            else if( l == 1 ) then
              cnd1 = .true.
              cnd2 = .true.
            end if
          end do

          if( cnd2 ) then
            c = 0.0
            s = 1.0
            do i = l,k
              f = s*rv1(i)
              if( (abs(f) + anorm) /= anorm ) then
                g = w(i)
                h = sqrt(f*f + g*g)
                w(i) = h
                h = 1.0/h
                c = (g*h)
                s = -(f*h)
                do j = 1,m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
                end do
              end if
            end do
          end if

          if( cnd1 ) then
            z = w(k)
            if( l == k ) then
              if( z < 0.0 ) then
                w(k) = -z
                do j = 1,n
                  v(j,k) = -v(j,k)
                end do
              end if
c             exit loop2
              go to 20
            end if
          end if

          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y)
          g = sqrt( f*f + 1.0 )
          f = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x
          c = 1.0
          s = 1.0
          do j = l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = sqrt( f*f + h*h )
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do nm = 1,n
              x = v(nm,j)
              z = v(nm,i)
              v(nm,j) = (x*c)+(z*s)
              v(nm,i) = -(x*s)+(z*c)
            end do
            z = sqrt( f*f + h*h )
            w(j) = z
            if( z /= 0.0 ) then
              z = 1.0/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do nm = 1,m
              y = a(nm,j)
              z = a(nm,i)
              a(nm,j) = (y*c)+(z*s)
              a(nm,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0
          rv1(k) = f
          w(k)   = x
        end do  !loop2
   20 continue
      end do
      
      end subroutine svdcmp

!-------------------------------------------------------------------------      
! purpose: solves a*x = b
!
! method:     
! solves a*x = b for a vector x, where a is specified by the arrays
! u,w,v as returned by svdcmp. m and n
! are the logical dimensions of a, and will be equal for square matrices.
! mp and np are the physical dimensions of a. b(1:m) is the input right-hand 
! side. x(1:n) is the output solution vector. no input quantities are 
! destroyed, so the routine may be called sequentially with different b
!
! author:  a. maute dec 2002   
! (* copyright (c) 1985 numerical recipes software -- svbksb *!
! from numerical recipes 1986 pp. 57 or can be find on web-sites
!-------------------------------------------------------------------------      

      subroutine svbksb( u, w, v, m, n, mp, np, b, x )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      integer, intent(in)   :: mp
      integer, intent(in)   :: np
      real , intent(in)  :: u(mp,np)
      real , intent(in)  :: w(np)
      real , intent(in)  :: v(np,np)
      real , intent(in)  :: b(mp)
      real , intent(out) :: x(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, j, jj
      real :: s
      real :: tmp(nmax)

      do j = 1,n
        s = 0. 
        if( w(j) /= 0. ) then
          do i = 1,m
            s = s + u(i,j)*b(i)
          end do
          s = s/w(j)
        endif
        tmp(j) = s
      end do

      do j = 1,n
        s = 0. 
        do jj = 1,n
          s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
      end do

      end subroutine svbksb
