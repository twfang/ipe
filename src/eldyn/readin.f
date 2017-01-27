!test022208: v2: how much degree can the polar cap conductanc affect the potential solution??? 
!
      module read_module
!      
      implicit none
      
      character(len=*),parameter :: input_type = 'ASCII'
!t      character(len=*),parameter :: input_type='NETCDF'

! working directory and filenames
!nm20140326
!#ifdef SUN
!      character(len=*),parameter :: work = 
!     &      '/suncat/e/maute/tiegcm_plasma_test/'	  ! (user given) 
!      character(len=*),parameter :: flnm  = 
!     &      'stest_org.nc' ! (user given)

!#elif AIX
!      character(len=*),parameter :: work = 
!     &      '/suncat/e/maute/tiegcm_plasma_test/'	  ! (user given) 
!      character(len=*),parameter :: flnm  = 
!     &      'stest_org.nc' ! (user given)
!#elif LINUX
      character(len=*),parameter :: work =                              &
     & '../tiegcm1.8_dynamo_hres/'  !nm082306: 
!nm082306:     & '/iapetus/i/naomi/tmp/tiegcm1.8_dynamo_hres/'  ! (user given) 
!nm040706:     &      '/suncat/e/maute/tiegcm_plasma_test/'	  ! (user given) 
      character(len=*),parameter :: flnm  =                             &
     &      'stest_NONE.nc' ! (user given)

!#endif
!
! ascii file
      character(len=*),parameter :: work_ascii =                        &
!nm20140403     &      '/suncat/e/maute/plasmasphere/'	  ! (user given) 
!nm20140408     &      './'
!nm20140730
     &''
!t     &'/home/Naomi.Maruyama/wamns/r311.1tmp/trunk/run/ipe_S_22711/'
!t     &'/home/Naomi.Maruyama/wamns/r319.2.3/trunk/run/ipe_S_30361/'
      character(len=*),parameter :: flnm_ascii  =                       &
!nm20140408     &      'integrals_plasma.inp' ! (user given)
!    &      'fort.4000' ! (user given)
     &      'stup_fli' ! (user given)
      character(len=*),parameter :: path_ascii=work_ascii//flnm_ascii
     
!nm20140403
      character(len=*),parameter :: work_ascii2 =                       &
     & '/home/Naomi.Maruyama/wamsc/bkup/naomilx/ctip_apex/dynamo/'

      !PRIVATE
      PUBLIC
!--------------------------------------------------------------------------  
      end module read_module
