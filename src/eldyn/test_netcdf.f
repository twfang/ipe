!module load intel netcdf 
!$ ifort -O3 -fp-model precise -r8 -I/apps/netcdf/3.6.3/intel/include test_netcdf.f -L/apps/netcdf/3.6.3/intel/lib -lnetcdf
      program test_netcdf
      implicit none
      INCLUDE 'netcdf.inc'
!        ...

      INTEGER  :: STATUS
      INTEGER  :: NCID
      INTEGER,parameter  :: NDIMS=3         ! number of dimensions
      INTEGER,parameter  :: TIMES=1  ! dimension lengths
      INTEGER,parameter  :: LATS=2   ! dimension lengths
      INTEGER,parameter  :: LONS=3  ! dimension lengths
      INTEGER  :: LATDIM,LONDIM,TIMDIM
      INTEGER  :: RHID               ! variable ID
      INTEGER  :: RHDIMS(3)               ! 
      INTEGER  :: START(NDIMS)
      INTEGER  :: COUNT(NDIMS)
      REAL ::    RHVALS(LONS, LATS, TIMES)
      DATA START /1, 1, 1/        ! start at first value
      DATA COUNT /LONS, LATS, TIMES/
      character(len=100) :: msgerr
      character :: fname*10,labl*56,units*12
      integer :: nw,nw_u
      integer :: ILON,ILAT,ITIME
!        ...

! create output file
      STATUS = nf_create('foo.nc',NF_CLOBBER,ncid)
      msgerr='nf_create'
      IF (STATUS .NE. NF_NOERR) CALL check_err(STATUS,msgerr)

 !     STATUS = NF_OPEN ('foo.nc', NF_WRITE, NCID)
 !     msgerr='nf_open'
 !     IF (STATUS .NE. NF_NOERR) CALL 
!!     & HANDLE_ERR(STATUS)
!     & check_err(STATUS)
!        ...

                                     ! define dimensions
      STATUS = NF_DEF_DIM(NCID, 'lat', 2, LATDIM)
      IF (STATUS .NE. NF_NOERR) CALL check_ERR(STATUS,msgerr)
      STATUS = NF_DEF_DIM(NCID, 'lon', 3, LONDIM)
      IF (STATUS .NE. NF_NOERR) CALL check_ERR(STATUS,msgerr)
      STATUS = NF_DEF_DIM(NCID, 'time', NF_UNLIMITED, TIMDIM)
      IF (STATUS .NE. NF_NOERR) CALL check_ERR(STATUS,msgerr)
 !             ...
                                      ! define variable


! end definition part      
      status = nf_enddef(ncid)
      if (status /= NF_NOERR) call check_err(status,'(1)enddef ')


      print *, 'ncid ', ncid 
       fname='rh'
       labl = 'rh'
       nw=2
       units = 'meter'
       nw_u=5

      STATUS = NF_INQ_VARID (NCID, fname, RHID)
      msgerr='NF_INQ_VARID'

!      IF (STATUS .NE. NF_NOERR) CALL 
!!HANDLE_ERR(STATUS)
!     & check_err(STATUS,msgerr)

      if (status /= NF_NOERR) then
        
        status = nf_redef(NCID)
        if (status /= NF_NOERR) call check_err(status,'redef ')

        print *,'(1) here after redef'
        RHdims(1) = LONDIM
        RHdims(2) = LATDIM
        RHdims(3) = TIMDIM
        status = nf_def_var(ncid,'rh',NF_DOUBLE,3,RHdimS,rhid)
        msgerr = 'defvar '//fname
        if (status /= NF_NOERR) call check_err(status,msgerr)
!
!      print *, '!dbg022007!a',noid,id
!
        status = nf_put_att_text(ncid,rhid,'long_name',nw,labl)
!      print *, '!dbg022007!b'

        msgerr = 'put att_'//fname//' name'

!        print *, '!dbg022007!c'

        if (status /= NF_NOERR) call check_err(status,msgerr)
!        print *, '!dbg022007!d',noid,id,nw_u,units      
        status = nf_put_att_text(ncid,rhid,'units',nw_u,units)
!        print *, '!dbg022007!e',istat
        msgerr = 'put att_'//fname//' units'
!        print *, '!dbg022007!f',istat,msgerr
        if (status /= NF_NOERR) call check_err(status,msgerr)
!        print *, '!dbg022007!g',noid
        status = nf_enddef(ncid)
!        print *, '!dbg022007!h'
        if (status /= NF_NOERR) call check_err(status,'(2)enddef ')
!        print *, '!dbg022007!i' 
      endif



      DO 10 ILON = 1, LONS
        DO 10 ILAT = 1, LATS
           DO 10 ITIME = 1, TIMES
              RHVALS(ILON, ILAT, ITIME) = 0.5
   10 CONTINUE
      print *,rhvals
      STATUS = NF_PUT_VARA_DOUBLE (NCID, RHID, START, COUNT, RHVALS)
      print *,' status', status,' NF_NOERR', NF_NOERR
      msgerr='NF_PUT_VARA_DOUBLE'
      IF (STATUS .NE. NF_NOERR) CALL 
!HANDLE_ERR(STATUS)
     & check_err(STATUS,msgerr)


      print *,'end program test_netcdf'
      end program test_netcdf

      subroutine check_err(status,name)
      implicit none
!
!nm20140326
!#ifdef SUN
!#include "/opt/local/include/netcdf.inc"
!#else
!#include "netcdf.inc"
      include "netcdf.inc"
!#endif
!
      integer,intent(in) :: status
      character(len=*),intent(in) :: name
!
      write(6,*) 'error in ',name
      write(6,*) 'NF_STRERROR=', NF_STRERROR(status)
!      write(6,*) ' '
!
      stop
! 
      end subroutine check_err

