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
! update integrals                                                                                                              
      module module_update_fli
      private
      public :: update_fli
      contains
      !---
      subroutine update_fli ( utime )
      use module_precision
      use module_input_parameters,ONLY:start_time,sw_debug,mype
      use dynamo_module,only:zigm11,zigm22,zigmc,zigm2,rim
      use module_plas2dyn_fli_array,only:plas2dyn_fli_array
      use params_module,only: kmlat,kmlonp1
      implicit none
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]      
!--local---
!NETCDF or IPE
      CHARACTER (len=*), parameter :: input_type_eld_fli = 'IPE'
!0: WITH wind; 1: WITHOUT wind
      INTEGER (KIND=int_prec) :: sw_windyn=1
      INTEGER (KIND=int_prec) :: i,j
!---
!note: when called at the very first time step, zigm11etc do not have values since IPE has not yet been called, therefore the readin values are used even if input_type_eld_fli='IPE'.

      if ( utime==start_time .or. input_type_eld_fli == 'NETCDF' ) then

!       if ( sw_debug )  
        print *, 'zigm11 from READIN: utime=',utime

! update every time step                                                                                                        
!        zigm11(:,:)=zigm11_readin(:,:)
!        zigm22(:,:)=zigm22_readin(:,:)
!        zigmc(:,:)=zigmc_readin(:,:)
!        zigm2(:,:)=zigm2_readin(:,:)
!        rim(:,:,:)=rim_readin(:,:,:)

      else if ( input_type_eld_fli == 'IPE' ) then

        print *, 'zigm11 from IPE: utime=',utime

! update values from IPE every time step                                                                                    
        call plas2dyn_fli_array ( utime )


!dbg20150608: output FLI for debug
        print *,'(22-26) output dyn fli at utime=',utime
!SMS$serial begin
        write(unit=4022,FMT='(I12)')utime
        write(unit=4022,FMT='(20E12.4)')zigm11
        write(unit=4023,FMT='(20E12.4)')zigm22
        write(unit=4024,FMT='(20E12.4)')zigmc
        write(unit=4025,FMT='(20E12.4)')zigm2
        write(unit=4026,FMT='(20E12.4)')rim
!SMS$serial end

      else
        write(6,*)'STOP! INVALID input_type_eld_fli=',input_type_eld_fli
        write(6,*)'mype=',mype
        stop
      endif    !( input_type_eld_fli == 'NETCDF' ) then                                                                           

! test with/without wind effect                                                                                                 
      if ( sw_windyn==0 ) then
        print *, 'update_fli: NO WIND'
        rim(:,:,:) = 0.0
      endif

!     if ( .NOT. sw_output_netcdf_eldyn ) then
!       print *,'output binary2: zigm11 etc: unit='
!     &     ,LUN_eldyn
!      write( UNIT = LUN_eldyn, IOSTAT=istat_eldyn ) zigm11
!     &     ,zigm22,zigm2,zigmc,rim
!      IF (istat_eldyn /= 0) STOP 'ERROR WRITING FILE_eldyn'
!      endif                     !( .NOT. sw_output_netcdf_eldyn ) then                                                          

!SMS$IGNORE BEGIN
      print *,'!dbg! make sure zigm11 have values before calculations'
      print *,'zigm11',mype,MAXVAL(zigm11),MINVAL(zigm11)
      print *,'zigm22',mype,MAXVAL(zigm22),MINVAL(zigm22)
      print *,'zigm2' ,mype,MAXVAL(zigm2 ),MINVAL(zigm2 )
      print *,'zigmc' ,mype,MAXVAL(zigmc ),MINVAL(zigmc )
      print *,'rim'   ,mype,MAXVAL(rim   ),MINVAL(rim   )
!SMS$IGNORE END

      do j=1,kmlat
        do i=1,kmlonp1
          if( zigm11(i,j) <= 0.0 ) then
!SMS$IGNORE begin
            print *,'!STOP! Bad zigm11 value in module_update_fli.f'
            print *,'VALUE',mype,i,j,zigm11(i,j)
!SMS$IGNORE end
            STOP
          endif
        enddo
      enddo

      end subroutine update_fli
      end module module_update_fli
