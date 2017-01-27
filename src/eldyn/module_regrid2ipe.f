!nm20140409 originally copied from impose2/readin.F
!nm101306: originally copied from /tiegcm1.8_dynamo_hres/
!         ncplot --> ncplot3D
!
      module module_regrid2ipe   
        implicit none
        PRIVATE 
        PUBLIC :: regrid2ipe         
      contains
!-----------------------------------------------------------------------
        subroutine regrid2ipe(ed_dyn,ed_IPE)
          USE module_precision
          USE module_IPE_dimension,ONLY: NPTS2D,NMP,NLP
          USE module_FIELD_LINE_GRID_MKS,ONLY:plasma_grid_GL,           &
     &        JMIN_IN,JMAX_IS
          use params_module,only:kmlon,kmlonp1,kmlat
          USE module_physical_constants,ONLY:rtd
          use cons_module,only: xlatm_deg
!---INTENT(IN)
!     INTEGER, DIMENSION(NLP),INTENT(IN) :: inGIP1d, isGIP1d ! 1d in, read from the IPE grid file     
!     REAL (KIND=8), DIMENSION(NPTS2D),INTENT(IN) :: ipeGL ! read from IPE grid file, but not used
      real(KIND=real_prec),INTENT(IN)::ed_dyn(kmlonp1,kmlat,2)
      real(KIND=real_prec),DIMENSION(2,NLP,NMP,2),INTENT(OUT) :: ed_IPE !1:N/SH....4:ed1/2
! Local 
      INTEGER (KIND=int_prec) :: IN,IS
      integer :: jth,jh,jl,ihem,imlat_ipe,imlat_dyn,imlat_dyn1    
      real :: fac,mlat_ipe
      real,dimension(2,kmlonp1,NLP,2) :: ed_lat ! interpolated in latitude for ipe

!     
! 
! mapping low resolution to high resolution
! interpolate in latitude     
      do jh = 1,NLP  ! loop over high resolution IPE
         IN = JMIN_IN(jh)
         IS = JMAX_IS(jh)
    
        do ihem=1,2
           if (ihem==1) then !NH
!             mlat_ipe = ipeGL(inGIP1d(jh))
              mlat_ipe = 90.-plasma_grid_GL(IN,jh)*rtd
           else if ( ihem==2 ) then !SH
!             mlat_ipe = ipeGL(isGIP1d(jh))
              mlat_ipe = 90.-plasma_grid_GL(IS,jh)*rtd
           end if
           imlat_ipe = INT(mlat_ipe*10.)

!        print *,'debug',ihem,IN,IS,mlat_ipe

        do jl = 1,kmlat  ! loop over low resolution:dyn
          imlat_dyn  = INT(xlatm_deg(jl  )*10.)
          imlat_dyn1 = INT(xlatm_deg(jl+1)*10.)
!	
	   if( imlat_ipe == imlat_dyn ) then
            do jth=1,2
	      ed_lat(ihem,:,jh,jth) = ed_dyn(:,jl,jth)
            end do 
!          print *,'chosen',jh,jl,imlat_ipe,imlat_dyn

	      goto 22
!	      
	   else if( imlat_ipe > imlat_dyn.and.imlat_ipe < imlat_dyn1 ) then
	     if(jl == kmlat) then
	       write(6,*) 'readin_netcdf: jl > kmlat',jl,kmlat
	     endif
	     fac = (mlat_ipe-xlatm_deg(jl))/                                   &
     &            (xlatm_deg(jl+1)-xlatm_deg(jl))
             do jth=1,2
	      ed_lat(ihem,:,jh,jth)= ed_dyn(:,jl,jth) + (ed_dyn(:,jl+1,jth)    &
     &          -ed_dyn(:,jl,jth))*fac
             end do 

	      
	      goto 22
	   endif !	   if( mlat_ipe == xlatm_deg(jl)) then
!	   
        end do !        do jl = 1,kmlat  ! loop over low resolution
  22	continue
        end do !ihem=1,2
      end do !      do jh = 1,NLP  ! loop over high resolution IPE
      
! interpolate in longitude is not needed, but rather need shifting ...
      do jh = 1,nmp  ! loop over high resolution IPE
         jl = jh + (nmp/2) !dynamo low resolution grid
         if (jl>nmp+1) jl=jl-nmp
!	
         do ihem=1,2
            do jth=1,2
               ed_IPE(ihem,:,jh,jth) = ed_lat(ihem,jl,:,jth)
            end do !jth
         end do                    !ihem
     
!	     
      enddo !      do jh = 1,nmp  ! loop over high resolution IPE 
      
      end subroutine regrid2ipe
      end module module_regrid2ipe
