!20111107: copied originally from tiegcm1.8_dynamo_lres
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      module module_sub_dynamo
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: dynamo
      contains
!-----------------------------------------------------------------------
!
! BOP
! !IROUTINE: dynamo 
! !INTERFACE:
      subroutine dynamo
!     !USES:
      use params_module,only: 
     |  kmlon,  ! number of geomagnetic grid longitudes
     |  kmlonp1,! kmlon+1
     |  kmlat,  ! number of geomagnetic grid latitudes
     |  kmlatp1,! kmlat+1
     |  kmlath  ! (kmlat+1)/2 (index to magnetic equator)
      use dynamo_module,only:zigm11,zigm22,zigmc,zigm2,rim,phim,phihm   &
     &,isolve,ncee,cee,nc,nc0,nc1,nc2,nc3,nc4,cofum,rhs,c0,c1,c2,c3,c4  &
     &,kmlon0,kmlon1,kmlon2,kmlon3,kmlon4,kmlat0,kmlat1,kmlat2,kmlat3   &
     &,kmlat4,pfrac,jn,jp
      use cons_module,only: dlatm,dlonm,pi_dyn,xlatm,rtd,jsecs      
      use nc_module,only: noid,	                                        &! id of output netcdf-file
     &     start3_out,                                                  &! only for put out 3D fields 
     &     dim3,count3
      use module_transf,ONLY:transf
      use module_rhspde,ONLY:rhspde
      use module_clearcee,ONLY:clearcee
      use module_stencmd,ONLY:stencmd
      use module_stencil,ONLY:stencil
      use module_edges,ONLY:edges
      use module_divide,ONLY:divide
      use module_stenmd,ONLY:stenmd
      use module_stenmod,ONLY:stenmod
!t      use module_mud,ONLY:mud
!t      use module_muh,ONLY:muh
!t      use module_mudmod,ONLY:mudmod
      use module_threed,ONLY:threed
!t      use module_sub_ncplot,ONLY:ncplot
      implicit none
!
! !DESCRIPTION:
! Transform needed fields to geomagnetic coordinates
! Perform field-line integrations
! Evaluate PDE coefficients and RHS
! The PDE is divided by 1/ DT0DTS (in dyncal divided by 1/cos(theta_0)
! Sigma_(phi phi) = zigm11/ rcos0s * dt0dts
! Sigma_(lam lam) = zigm22 * rcos0s / dt0dts
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     =   rim(1) * dt0dts
! K_(m lam)^D     = +-rim(2) * rcos0s
!
! !RETURN VALUE:
! !PARAMETERS: 
! !REVISION HISTORY:
! 05.02.04  <Astrid Maute> <include header> 
! 
! EOP
! 
! Local:
      integer :: i,j,jj,jjj,j0,jntl,k,n,ncc,nmaglat,nmaglon,ier
      real :: sym
      real :: array(-15:kmlon0+16,kmlat0),cs(kmlat0)
      character :: fname*10,labl*56,units*12
     

      print *, 'sub-dyn: isolve=' ,isolve 
!

      call transf
!
! Fold southern hemisphere over on to northern  (was in transf.F version tgcm15)
! -Value at the equator is also folded therefore at the equatorial boundary
!  factor of 1/2 introduced
! -added values in array index 49-97 from equator to poles
! -SH array (1-49) has original SH values
! -reverse sign of zigmc (to be compatible with Cicely's) 
! -sign of K_(m lam)^D in southern hemisphere is reversed, therefore 
!     K_(m lam)^N - (- K_(m lam)^S)
!
! zigm11 = Sigma_(phi phi)(0)^T
! zigm22 = Sigma_(lam lam)(0)^T
! zigmc  =-Sigma_c(0)^T 
! zigm2  = Sigma_h(0)^T 
! rim(1) = K_(m phi)^D(0)^T
! rim(2) = K_(m lam)^D(0)^T
!




      do j=1,kmlath
        do i=1,kmlonp1
          zigm11(i,kmlatp1-j) = (zigm11(i,kmlatp1-j)+zigm11(i,j))
          zigmc(i,kmlatp1-j)  =  zigmc(i,kmlatp1-j)+zigmc(i,j)
          zigmc(i,kmlat+1-j)  = -zigmc(i,kmlat+1-j)
          zigm2(i,kmlat+1-j)  =  zigm2(i,kmlat+1-J)+zigm2(i,j)
          zigm22(i,kmlat+1-j) = (zigm22(i,kmlat+1-J)+zigm22(i,j))
	  
          rim(i,kmlatp1-j,1)  = (rim(i,kmlat+1-j,1)+rim(i,j,1))
          rim(i,kmlatp1-j,2)  = (rim(i,kmlat+1-j,2)+rim(i,j,2))  
        enddo ! i=1,kmlonp1
      enddo ! j=1,kmlath
!
! Tzu-Wei TEST
      print *,'TEST zigm within dynamo'
      print *,'zigm11',MAXVAL(zigm11),MINVAL(zigm11)
      print *,'zigm22',MAXVAL(zigm22),MINVAL(zigm22)
      print *,'zigm2',MAXVAL(zigm2),MINVAL(zigm2)
      print *,'zigmc',MAXVAL(zigmc),MINVAL(zigmc)
      print *,'rim',MAXVAL(rim),MINVAL(rim)

! Calculate RHS of PDE from rim(1) and rim(2) (was in transf.F version tgcm15)
!  [( d K_(m phi)^D / d phi /(cos(theta_m)?) +
!  (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) ] * R / (RCOS0S*DT0DTS)
! ~ J_(Mr)*r^2*cos(theta_m)/cos(theta_0)/DT0DTS
!

      call rhspde

!
! Set index array nc and magnetic latitude cosine array:
! nc pointes to the start of the coefficient array for each level
      nc(1) = nc0
      nc(2) = nc1
      nc(3) = nc2
      nc(4) = nc3
      nc(5) = nc4
      nc(6) = ncee
!
! Use pi_dyn from cons module rather than 4*atan(1) to avoid small
! differences generated by the atan in -lmass lib (-lmass was not
! used in earlier versions).
!
      do j=1,kmlat0
        cs(j) = cos(pi_dyn/2.-(kmlat0-j)*dlatm)
      enddo ! j=1,kmlat0
!
! Set up difference coefficients. 
! mixed terms: factor 4 from 5-point diff. stencil
! "NH" array index 49:97 equator to pole
! zigmc  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
!
! "SH" array index 49:1 equator to pole
! zigmc  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
! 
!
      j0 = kmlat0-kmlath
      do j=1,kmlath       !  1,49 (assuming kmlat=97)
        jj = kmlath+j-1   ! 49,97 added values ()^T 
        jjj = kmlath-j+1  ! 49,1
!
        do i=1,kmlonp1
          zigmc(i,jj)   = (zigmc(i,jj) +zigm2(i,jj))/(4.*dlatm*dlonm)
          zigm2(i,jj)   = zigmc(i,jj)-2.*zigm2(i,jj)/(4.*dlatm*dlonm)
          zigm22(i,jj)  = zigm22(i,jj)*cs(j0+j)/dlatm**2
          zigmc(i,jjj)  = -zigmc(i,jj)
          zigm2(i,jjj)  = -zigm2(i,jj)
          zigm22(i,jjj) = zigm22(i,jj)
        enddo ! i=1,kmlonp1
        if (j /= kmlath) then
          do i = 1,kmlonp1
            zigm11(i,jj) = zigm11(i,jj)/(cs(j0+j)*dlonm**2)
            zigm11(i,jjj) = zigm11(i,jj)
          enddo
        endif
      enddo ! j=1,kmlath


!
! Set zigm11 to zero at megnetic poles to avoid floating exception 
! (values at poles are not used):
!
      do i = 1,kmlonp1
        zigm11(i,1)     = 0.0
        zigm11(i,kmlat) = 0.0
      enddo
!
! Clear array for difference stencils at all levels:
      call clearcee(cee,kmlon0,kmlat0)

!
! Calculate contribution to stencils from each PDE coefficient
! isolve = 0 -> original mud version 5.
! isolve = 1 -> muh hybrid solver (only as direct solver -- slow)
! isolve = 2 -> modified mudpack solver (modified and unmodified coefficients)
!
      if (isolve==2) then
        cofum(:,:,:) = 0. ! init
!
! Sigma_(phi phi)(0)/( cos(lam_0)*(d lon)^2 )
        call stencmd(zigm11(1,kmlat0),kmlon0,kmlat0,cee,1)

!
! Sigma_(lam lam)(0)*cos(lam_0)/(d lam_0)^2
        call stencmd(zigm22(1,kmlat0),kmlon0,kmlat0,cee,4)

!
! Sigma_(phi lam)(0)/( 4*d lam_0* d lon )
	zigmc(:,kmlat0) = -zigmc(:,kmlat0)
        call stencmd(zigmc(1,kmlat0),kmlon0,kmlat0,cee,2)
	zigmc(:,kmlat0) = -zigmc(:,kmlat0)

!
! Sigma_(lam phi)(0)/( 4*d lam_0* d lon )
	zigm2(:,kmlat0) = -zigm2(:,kmlat0)
        call stencmd(zigm2(1,kmlat0),kmlon0,kmlat0,cee,3)
	zigm2(:,kmlat0) = -zigm2(:,kmlat0)

!
! isolve /= 2: original or hybrid solver (only modified stencil).
      else
!
! Sigma_(phi phi)(0)/( cos(lam_0)*(d lon)^2 )
        call stencil(zigm11(1,kmlat0),kmlon0,kmlat0,cee,1)

! Sigma_(lam lam)(0)*cos(lam_0)*/(d lam_0)^2
        call stencil(zigm22(1,kmlat0),kmlon0,kmlat0,cee,4)

!
! Sigma_(phi lam)(0)/( 4*d lam_0* d lon )
	zigmc(:,kmlat0) = -zigmc(:,kmlat0)
        call stencil(zigmc(1,kmlat0),kmlon0,kmlat0,cee,2)
	zigmc(:,kmlat0) = -zigmc(:,kmlat0)

!
! Sigma_(lam phi)(0)/( 4*d lam_0* d lon )
	zigm2(:,kmlat0) = -zigm2(:,kmlat0)
        call stencil(zigm2(1,kmlat0),kmlon0,kmlat0,cee,3)
	zigm2(:,kmlat0) = -zigm2(:,kmlat0)

!
      endif ! isolve
!
! Insert RHS in finest stencil (formerly sub rths):
      do j = 1,kmlat0
        jj = kmlath-kmlat0+j
        do i = 1,kmlon0
          c0(i,j,10) = rhs(i,jj)
        enddo ! i = 1,kmlon0
      enddo ! j = 1,kmlat0
!
! Set boundary condition at the pole:
      call edges(c0,kmlon0,kmlat0)
      call edges(c1,kmlon1,kmlat1)
      call edges(c2,kmlon2,kmlat2)
      call edges(c3,kmlon3,kmlat3)
      call edges(c4,kmlon4,kmlat4)
      if (isolve==2) 
     |  call edges(cofum,kmlon0,kmlat0)
!
! Divide stencils by cos(lam_0) (not rhs):
      call divide(c0,kmlon0,kmlat0,kmlon0,kmlat0,cs,1)
      call divide(c1,kmlon1,kmlat1,kmlon0,kmlat0,cs,1)
      call divide(c2,kmlon2,kmlat2,kmlon0,kmlat0,cs,1)
      call divide(c3,kmlon3,kmlat3,kmlon0,kmlat0,cs,1)
      call divide(c4,kmlon4,kmlat4,kmlon0,kmlat0,cs,1)
      if (isolve==2) 
     |  call divide(cofum,kmlon0,kmlat0,kmlon0,kmlat0,cs,0)
!
! Set value of solution to 1. at pole:
      do i=1,kmlon0
        c0(i,kmlat0,10) = 1.
      enddo
!
! Modify stencils and RHS so that the NH high lat potential is inserted at
!  high latitude.  The SH high lat potential will be added back later.
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
!
      ncc = 1
      nmaglon = kmlon0
      nmaglat = kmlat0
      do n=1,5
        if (isolve==2) then
          call stenmd(nmaglon,nmaglat,cee(ncc),phihm(1,kmlat0),pfrac)
        else
          call stenmod(nmaglon,nmaglat,cee(ncc),phihm(1,kmlat0),pfrac)
        endif
        ncc = ncc+9*nmaglon*nmaglat
        if (n==1) ncc = ncc+nmaglon*nmaglat ! rhs is in 10th slot
        nmaglon = (nmaglon+1)/2
        nmaglat = (nmaglat+1)/2
      enddo ! n=1,5
!
      jntl = 0
!
      ier = 0 
      if(isolve==0) then
        call mud(rim,jntl,isolve,ier)	 ! solver in mud.F
        if(ier < 0 ) then 	! not converged
	  write(6,*) 'muh: use direct solver'
	  call muh(rim,jntl) 		! solver in mud.F
	endif
      elseif (isolve==1) then
        call muh(rim,jntl)        	! solver in muh2cr.F
      elseif (isolve==2) then
        call mudmod(rim,jntl,isolve,ier)! solver in mudmod.F
        if(ier < 0 ) then 	! not converged
	  write(6,*) 'muh: use direct solver'
	  call muh(rim,jntl) 		! solver in mud.F
	endif
      else
        write(6,*) 'dynamo: solver type ',isolve,' not implemented.'
        stop 'isolve'
      endif
!
! Copy output potential from rim to phim(kmlonp1,kmlat):
!  Correct the SH potential for the anti-symmetric imposed NH high lat poten
!  pfrac = fraction of dynamo in solution in the NH. = 1 low lat, = 0 hi lat
!    cons_module: crit(1)=15, crit(2)=30 deg colats, or hi-lat > 75 deg,
!      dynamo < 60 deg, and combination between 60-75 mag lat.
!  ie, need (1.-pfrac(i,kmlat0)) at phim(i,1), etc
! jn index for NH part of potential (kmlat down to ~kmlat0)
! jp index for NH pfrac (kmlat0 down to 1)
!
      do j=1,kmlat0
        jn = kmlat - j + 1
        jp = kmlat0 - j + 1
        do i=1,kmlonp1
          phim(i,j)=rim(i,j,1)+(1.-pfrac(i,jp))*(phihm(i,j)-phihm(i,jn))
!	  write(6,'(2i4,3(x,f12.5))') i,j,rim(i,j,1),pfrac(i,jp),
!     |          phihm(i,j)
        enddo ! i=1,kmlonp1
      enddo ! j=1,kmlat0
      do j=kmlat0+1,kmlat
        jp = j - kmlat0
        do i=1,kmlonp1
          phim(i,j) = rim(i,j,1)
        enddo ! i=1,kmlonp1
      enddo ! j=1,kmlat 

      fname = 'poten'
      labl = 'poten'
      units = 'V'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,
     |  	  phim,5,units,1)

!dbg20140801(27)
      print *,'(27)output dyn PHIM at utime=',jsecs
!      write(unit=4027,FMT='(I12)')INT(jsecs)
      write(unit=4027,FMT='(20E12.4)')phim

!
! Call threed to calculate 2-d electric potential array in geomagnetic coordinates
!   from 2-d solver output phim, corrected for the SH potential
!
      call threed
!    
      end subroutine dynamo
!-----------------------------------------------------------------------
      end module module_sub_dynamo
!-----------------------------------------------------------------------
