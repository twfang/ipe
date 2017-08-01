!test022208: v2: how much degree can the polar cap conductanc affect the potential solution??? 
!
      module module_readin_ascii      
!      
      PRIVATE
      PUBLIC :: readin_ascii
      contains
!-----------------------------------------------------------------------
      subroutine readin_ascii
      use params_module,only:                                           &
     &  kmlonp1,                                                        &! kmlon+1
     &  kmlat   ! number of geomagnetic grid latitudes

      use read_module,only:work_ascii2,path_ascii

! read in data file with integrals from George Millward
!
      use dynamo_module,only:                                           &
     &  zigm11,                                                         &! zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
     &  zigm22,                                                         &! zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
     &  zigmc,                                                          &! zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
     &  zigm2,                                                          &! zigm2  is int[sigma_h] ds,	      i.e. Sigma_h
     &  rim      ! K_(m phi)^D/abs(sin Im),  K_(m lam)^D
      
      use nc_module,only: noid,	                                        &! id of output netcdf-file
     &     start3_out,                                                  &! only for put out 3D fields 
     &     dim3,count3
!t      use module_sub_ncplot,ONLY:ncplot
!t      use cons_module,only: secs
      implicit none
!
! Local 
      integer:: lp,mp,i,is,ie,itotal,j
      integer, parameter :: iunit= 179,ounit=4001
      character :: fname*10,labl*56,units*12
      real,dimension(kmlonp1,kmlat) ::                                  &
     &  part1,part2,part3  ! for comparison with k_mlam terms
      integer,parameter:: sw_part1=0
      integer :: jth,utime
      real,dimension(kmlonp1,kmlat)::dum
!      
!SMS$SERIAL BEGIN
      open(iunit,file=path_ascii,status='OLD',err=99)
      write(6,*) 'open ascii file ',path_ascii
      read(unit=iunit,FMT="(I12)",err=199)utime
      print *,'utime=',utime
!t      secs = REAL(utime)
!t      print *,'secs=',secs
!
      jth_loop: do jth=1,6
      write(6,*) 'read in jth=',jth
        read(unit=iunit,FMT="(20E12.4)")dum





        if (jth==1) then 
           zigm11=dum
        else if (jth==2) then 
           zigm22=dum
        else if (jth==3) then 
           zigm2 =dum
        else if (jth==4) then 
           zigmc =dum
        else if (jth==5) then 
           rim(:,:,1)=dum
        else if (jth==6) then 
           rim(:,:,2)=dum
        end if
      end do jth_loop


!dbg20140804: bug is fixed in IPE so i dont need it any more
!      if ( MINVAL(zigm11) <= 0.0 ) then 
!         print *,'!STOP! read ascii: INVALID MIN zigm11', MINVAL(zigm11)
!!         STOP
!         zigm11( 1,:)=( zigm11(2,:)+zigm11(80,:) )*0.5
!         zigm11(81,:)=zigm11(1,:)
!         zigm22( 1,:)=( zigm22(2,:)+zigm22(80,:) )*0.5
!         zigm22(81,:)=zigm22(1,:)
!         zigm2( 1,:)=( zigm2(2,:)+zigm2(80,:) )*0.5
!         zigm2(81,:)=zigm2(1,:)
!         zigmc( 1,:)=( zigmc(2,:)+zigmc(80,:) )*0.5
!         zigmc(81,:)=zigmc(1,:)
!         do jth=1,2
!            rim( 1,:,jth)=( rim(2,:,jth)+rim(80,:,jth) )*0.5
!            rim(81,:,jth)=rim(1,:,jth)
!         end do
!      endif

!
!nm20140408--start commented out
!      itotal = 0
!      write(6,*) 'read in zigm11'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!          read(iunit,2435) (zigm11(mp,lp),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (zigm11(mp,lp))
      
!      write(6,*) 'read in zigm22'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!     read(iunit,2435) (zigm22(mp,lp),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (zigm22(mp,lp))
!      enddo
      
!      write(6,*) 'read in zigm2'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!          read(iunit,2435) (zigm2(mp,lp),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (zigm2(mp,lp))
!      enddo
      
!      write(6,*) 'read in zigmc'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!          read(iunit,2435) (zigmc(mp,lp),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (zigmc(mp,lp))
!      enddo
      
!      write(6,*) 'read in rim1'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!          read(iunit,2435) (rim(mp,lp,1),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (rim(mp,lp,1))
!      enddo
      
!      write(6,*) 'read in rim2'
!      do lp = 1,kmlat
!        is = -9
!        do i=1,8
!	  is = is + 10
!	  ie = is+9
!          read(iunit,2435) (rim(mp,lp,2),mp=is,ie)
!	enddo
!	mp = ie + 1
!	itotal = itotal + mp
!        read(iunit,2436) (rim(mp,lp,2))
!      enddo
!      write(6,*) 'itotal ',itotal
!
!nm20140408--end commented out
      close(iunit,err=100)
!   
! for testing input: is the same as input file
!
!      open(ounit,file='/suncat/e/maute/plasmasphere/read.out',err=99)
!      do lp = 1,kmlat
!	  write(ounit,2435) (zigm11(mp,lp),mp=1,kmlonp1)
!      enddo
!      do lp = 1,kmlat
!          write(ounit,2435) (zigm22(mp,lp),mp=1,kmlonp1)
!      enddo
!      do lp = 1,kmlat
!          write(ounit,2435) (zigm2(mp,lp),mp=1,kmlonp1)
!      enddo
!      do lp = 1,kmlat
!          write(ounit,2435) (zigmc(mp,lp),mp=1,kmlonp1)
!      enddo
!      do lp = 1,kmlat
!          write(ounit,2435) (rim(mp,lp,1),mp=1,kmlonp1)
!      enddo
!      do lp = 1,kmlat
!          write(ounit,2435) (rim(mp,lp,2),mp=1,kmlonp1)
!      enddo
!      close(ounit,err=100)

	
! am 10/04 so far no value at the equator therefore set it
! but this shouldn't be in the code
      j = kmlat/2+1
      do i = 1,kmlonp1
	zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
	zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
	zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
	zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
	rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
	rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
      enddo ! i = 1,kmlon
      
! for testing
!      zigmc = 0.
!      zigm2 = 0.
!      rim(:,:,1) = 0
! am 10/04 testing positive K_mlam during night time
!        rim(1:30,:,2) = 0.
!        rim(70:81,:,2) = 0.

! am 10/04 change sign of K_(m lam)^D in the SH- that's what TIEGCM dynamo
! expects
      do lp = 1,(kmlat+1)/2
!        rim(:,lp,2) = rim(:,lp,2)*1.e8
        rim(:,lp,2) = -rim(:,lp,2)
      enddo
!      
! am 10/04 set high latitude values since these are no calculated
! in the plasmasphere model region of Heelis pattern
! value itself should not matter since potential is prescribed  
!t      do lp = 1,14
!t        zigm11(:,lp) = 0.01
!t        zigm22(:,lp) = 0.01
!t        zigm11(:,kmlat+1-lp) = 0.01
!t        zigm22(:,kmlat+1-lp) = 0.01
!t      enddo
!SMS$SERIAL END
      
! ouput variables     
      fname = 'zigm11_in'
      labl = 'zigm11_in'
      units = 'S'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  zigm11,9,units,1)  

!dbg20140801(1)
      print *,'output dyn fli at utime=',utime
      write(unit=ounit,FMT='(I12)')utime
      write(unit=ounit,FMT='(20E12.4)')zigm11

      fname = 'zigm22_in'
      labl = 'zigm22_in'
      units = 'S'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  zigm22,9,units,1)    
      fname = 'zigm2_in'
      labl = 'zigm2_in'
      units = 'S'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  zigm2,8,units,1)    
      fname = 'zigmc_in'
      labl = 'zigmc_in'
      units = 'S'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  zigmc,8,units,1)   
      fname = 'rim1_in'
      labl = 'rim1_in'
      units = 'A/m'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  rim(:,:,1),7,units,3) 
      fname = 'rim2_in'
      labl = 'rim2_in'
      units = 'A/m'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  rim(:,:,2),7,units,3)
     
!  
! test for k_mlam part1, part2,part3    
      if (sw_part1==1) then
      open(iunit,file=                                                  &
     & work_ascii2//'fort.179_4',                                       &!nm20140403:
!nm20140403     &'/suncat/e/maute/plasmasphere/part123.inp',
     &   status='OLD',err=100)
!      write(6,*) 'open file ',path
!     
      itotal = 0
      write(6,*) 'read in part1'
      do lp = 1,kmlat
        is = -9
        do i=1,8
	  is = is + 10
	  ie = is+9
          read(iunit,2435) (part1(mp,lp),mp=is,ie)
	enddo
	mp = ie + 1
	itotal = itotal + mp
        read(iunit,2436) part1(mp,lp)
      enddo
      write(6,*) 'read in part2'
      do lp = 1,kmlat
        is = -9
        do i=1,8
	  is = is + 10
	  ie = is+9
          read(iunit,2435) (part2(mp,lp),mp=is,ie)
	enddo
	mp = ie + 1
	itotal = itotal + mp
        read(iunit,2436) part2(mp,lp)
      enddo
      write(6,*) 'read in part3'
      do lp = 1,kmlat
        is = -9
        do i=1,8
	  is = is + 10
	  ie = is+9
          read(iunit,2435) (part3(mp,lp),mp=is,ie)
	enddo
	mp = ie + 1
	itotal = itotal + mp
        read(iunit,2436) part3(mp,lp)
      enddo
!      
      close(iunit,err=100)  
!      
      fname = 'part1_in'
      labl = 'part1_in'
      units = 'A/m'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  part1(:,:),8,units,3)
      fname = 'part2_in'
      labl = 'part2_in'
      units = 'A/m'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  part2(:,:),8,units,3)
      fname = 'part3_in'
      labl = 'part3_in'
      units = 'A/m'
      call ncplot(noid,fname,labl,start3_out,count3,dim3,               &
     &  	  part3(:,:),8,units,3)
      end if !(sw_part1==1) then
!      

!dbg20140801(2)
      print *,'output dyn fli at utime=',utime
      write(unit=4002,FMT='(I12)')utime
      write(unit=4002,FMT='(20E12.4)')zigm11

      return
!      
 2435 format(10e12.4)
 2436 format(e12.4)
   99 write(6,*) 'Error opening input file',path_ascii
      STOP
  199 write(6,*) 'Error reading input file'
      STOP
  100 write(6,*) 'Error closing input file'
!      
      end subroutine readin_ascii 
!--------------------------------------------------------------------------  
      end module module_readin_ascii
