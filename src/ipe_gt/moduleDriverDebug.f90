MODULE moduleDriverDebug

IMPLICIT NONE

PRIVATE

PUBLIC :: checkGridIPE
PUBLIC :: check2GridIPE
PUBLIC :: checkThermo
PUBLIC :: checkThermoArrays
PUBLIC :: checkInterp
PUBLIC :: checkFixedGeo

CONTAINS

!==========================================================================

SUBROUTINE checkGridIPE(directory, fileName, unitNumber, NMP, NLP, &
                        inGIP1d, isGIP1d, mlon_rad)

IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER :: unitNumber

INTEGER, intent(in) :: NLP, NMP
INTEGER, intent(in) :: inGIP1d(NLP)
INTEGER, intent(in) :: isGIP1d(NLP)
REAL (KIND=8), intent(in) :: mlon_rad(1:NMP+1)


!--------------------
! local variables 
!--------------------
!------------------------------------------------------
! File unit number for checking IPE grid values :
!------------------------------------------------------
INTEGER ::  FileOpenStatus  ! for checking if we can open the file
INTEGER :: ii

CHARACTER(100) :: fullFileName

! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!================================================
! Open for checking values from the IPE grid file
!================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkGridIPE : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

print *,'checkGridIPE : writing to ',fullFileName


WRITE (unitNumber, '(A)' ) 'inGIP1d : '
do ii = 1, NLP, 10
   WRITE (unitNumber, '(10I8)' ) inGIP1d(ii:ii+9)
enddo 

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A)' ) 'isGIP1d : '
do ii = 1, NLP, 10
   WRITE (unitNumber, '(10I8)' ) isGIP1d(ii:ii+9)
enddo 

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A)' ) 'mlon_rad : '
do ii = 1, NMP, 10
   WRITE (unitNumber, '(10ES15.5)' ) mlon_rad(ii:ii+9)
enddo 
WRITE (unitNumber, '(ES15.5)' ) mlon_rad(NMP+1)



WRITE (unitNumber, '(A,F15.5)') 'mlon_rad(2) = ', mlon_rad(2)
WRITE (unitNumber, '(A)') ' '

END SUBROUTINE checkGridIPE



! =======================================================================================


SUBROUTINE check2GridIPE(unitNumber, NPTS, NMP, &
                         glat_plasma_3d, glond_plasma_3d, pz_plasma_3d)

IMPLICIT NONE

INTEGER :: unitNumber

INTEGER, intent(in) :: NPTS, NMP


REAL(kind=8), intent(in) :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp)
REAL(kind=8), intent(in) :: pz_plasma_3d(npts,nmp)


!--------------------
! local variables 
!--------------------

INTEGER :: ii, jj 

! BEGIN CODE ====================================================================

ii = 1118
jj = 1

WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'Pz_plasma_3d(', ii, ',' , jj, ') = ', Pz_plasma_3d(ii,jj)

jj = 2
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLAT_plasma_3d(' , ii, ',' , jj, ') = ', GLAT_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLONd_plasma_3d(', ii, ',' , jj, ' ) = ', GLONd_plasma_3d(ii,jj)

WRITE (unitNumber, '(A)') ' '

ii = 700
jj =  20
!WRITE (unitNumber, '(A,F15.5)') 'Pz_plasma_3d(700,20) = ', Pz_plasma_3d(ii,jj)
!WRITE (unitNumber, '(A,F15.5)') 'GLAT_plasma_3d(700,20) = ', GLAT_plasma_3d(ii,jj)
!WRITE (unitNumber, '(A,F15.5)') 'GLONd_plasma_3d(700,20) = ', GLONd_plasma_3d(ii,jj)

WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'Pz_plasma_3d(' , ii, ',' , jj, ') = ', Pz_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLAT_plasma_3d(' , ii, ',' , jj, ') = ', GLAT_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLONd_plasma_3d(', ii, ',' , jj, ' ) = ', GLONd_plasma_3d(ii,jj)


END SUBROUTINE check2GridIPE


!=============================================================================================================

SUBROUTINE checkThermo(directory, fileName, unitNumber, &
       GT_input_dataset, GT_output_dataset, nday, &
       Universal_Time_seconds, solar_declination_angle_radians,  &
       houghSize, hough11, hough22, hough23, hough24, hough25, &
       amplitude, tidalPhase, GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT, Ne_density_FOR_GT)


USE moduleAmplitude  ! has amplitude type
USE moduleTidalPhase ! has tidalPhase type

IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber
CHARACTER(*), intent(in) :: GT_input_dataset , GT_output_dataset
INTEGER, intent(in) :: nday
REAL*8, intent(in) :: universal_time_seconds 
REAL*8, intent(in) :: solar_declination_angle_radians
INTEGER, intent(in) :: houghSize
REAL*8, intent(in) :: hough11(houghSize)
REAL*8, intent(in) :: hough22(houghSize)
REAL*8, intent(in) :: hough23(houghSize)
REAL*8, intent(in) :: hough24(houghSize)
REAL*8, intent(in) :: hough25(houghSize)
TYPE(amplitudeType), intent(in) :: amplitude
TYPE(tidalPhaseType), intent(in) :: tidalPhase

INTEGER, intent(in) :: GT_ht_dim, GT_lat_dim, GT_lon_dim
REAL(kind=8), intent(in) :: wind_southwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wind_eastwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wvz_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: rmt_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)

REAL*8, intent(in) :: Ne_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


! Local Variables ----------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file

! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!================================================
! Open for checking values from GT 
!================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkThermo : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

print *,'checkThermo : writing to ',fullFileName

WRITE (unitNumber, '(A)') GT_input_dataset
WRITE (unitNumber, '(A)') GT_output_dataset
WRITE (unitNumber, '(A,I5)') 'nday = ', nday
WRITE (unitNumber, '(A, F20.7)' ) 'Universal_Time_seconds = ',  Universal_Time_seconds
WRITE (unitNumber, '(A, 3x, F10.5)' ) 'solar_declination_angle_radians = ', &
                    solar_declination_angle_radians
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough11 = ',hough11
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough22 = ',hough22
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough23 = ',hough23
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough24 = ',hough24
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough25 = ',hough25

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 5F10.5)' ) 'amplitude = ', amplitude

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 5F10.5)' ) 'tidalPhase = ', tidalPhase


CALL checkThermoArrays(unitNumber, &
       GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT)


WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(1,1,:) = ', Ne_density_FOR_GT(1,1,:)
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(11,51,:) = ', Ne_density_FOR_GT(11,51,:)
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(15,91,:) = ', Ne_density_FOR_GT(15,91,:)


END SUBROUTINE checkThermo

!=================================================================================================



SUBROUTINE checkThermoArrays(unitNumber, &
       GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT)


IMPLICIT NONE

INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) :: GT_ht_dim, GT_lat_dim, GT_lon_dim
REAL(kind=8), intent(in) :: wind_southwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wind_eastwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wvz_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: rmt_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


! Local Variables ----------


! BEGIN CODE ====================================================================


WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(1,1,:) = ', &
                      wind_southwards_ms1_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(11,51,:) = ', &
                      wind_southwards_ms1_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(15,91,:) = ', &
                      wind_southwards_ms1_FROM_GT(15,91,:)
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(1,1,:) = ', &
                      wind_eastwards_ms1_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(11,51,:) = ', &
                      wind_eastwards_ms1_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(15,91,:) = ', &
                      wind_eastwards_ms1_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(1,1,:) = ', wvz_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(11,51,:) = ', wvz_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(15,91,:) = ', wvz_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(1,1,:) = ',rmt_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(11,51,:) = ',rmt_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(15,91,:) = ',rmt_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(1,1,:) = ',Temperature_K_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(11,51,:) = ',Temperature_K_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(15,91,:) = ',Temperature_K_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(1,1,:) = ',ht_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(11,51,:) = ',ht_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(15,91,:) = ',ht_FROM_GT(15,91,:)




END SUBROUTINE checkThermoArrays


!=================================================================================================


SUBROUTINE checkInterp(directory, fileName, unitNumber, &
                       interface_hts, GT_lat_dim, GT_lon_dim, &
                       O_density_fixed_ht, O2_density_fixed_ht, N2_density_fixed_ht, &
                       Vx_fixed_ht, Vy_fixed_ht, wvz_fixed_ht, &
                       tts_fixed_ht, qion3d_fixed_ht, elx_fixed_ht, &
                       ely_fixed_ht)


IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) ::  interface_hts, GT_lat_dim, GT_lon_dim


REAL(kind=8), intent(in) :: O_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: O2_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: N2_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Vx_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Vy_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Wvz_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: tts_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: qion3d_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: elx_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ely_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)



! Local Variables ---------------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file
INTEGER ::  ii, jj

! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!======================================================
! Open file for checking values from the interpolation
!======================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkInterp : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

!========================================
! Open for checking values from the interpolation
!========================================
!OPEN (UNIT=unitNumber, FILE='/mtb/save/wx20lrm/DATA/IONOSPHERE/gtgipIPE/OUTPUT/interpOut.dat', &
!      ACTION="WRITE", IOSTAT=FileOpenStatus)
!IF (FileOpenStatus > 0) STOP "Cannot open  interpOut.dat file ********"


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

       WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'O_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
       WRITE (unitNumber, '(20ES20.10)' ) O_density_fixed_ht(ii,jj,:)

    end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'O2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) O2_density_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'N2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) N2_density_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'Vx_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) Vx_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'Vy_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) Vy_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'wvz_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) wvz_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'tts_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) tts_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'qion3d_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) qion3d_fixed_ht(ii,jj,:)

    end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'elx_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) elx_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'ely_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) ely_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


END SUBROUTINE checkInterp

!====================================================================================================


SUBROUTINE checkFixedGeo(directory, fileName, unitNumber, &
                         npts, nmp, &
                         TN_plasma_input_3d, O_plasma_input_3d, O2_plasma_input_3d, &
                         N2_plasma_input_3d, GLAt_plasma_3d, GLOnd_plasma_3d, &
                         pz_plasma_3d, V_east_plasma, V_south_plasma, &
                         V_upward_plasma, ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
                         ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, ispecial_3d_fixed_ht, &
                         ihl_3d_fixed_ht, ihu_3d_fixed_ht, isFirstCallFixedHeight)


IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) :: npts, nmp

REAL(kind=8), intent(in) :: TN_plasma_input_3d(npts, nmp), O_plasma_input_3d(npts, nmp), &
                            O2_plasma_input_3d(npts, nmp), N2_plasma_input_3d(npts, nmp)

REAL(kind=8), intent(in)  :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp), pz_plasma_3d(npts,nmp)

REAL(kind=8), intent(in) :: V_east_plasma(npts, nmp), &
                            V_south_plasma(npts, nmp), &
                            V_upward_plasma(npts, nmp)

INTEGER, intent(in) :: ilon1_3d_fixed_ht(npts,nmp), ilon2_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ilat1_3d_fixed_ht(npts,nmp), ilat2_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ispecial_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ihl_3d_fixed_ht(npts,nmp), ihu_3d_fixed_ht(npts,nmp)
LOGICAL, intent(in)  :: isFirstCallFixedHeight

! Local Variables ---------------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file


! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!======================================================
! Open file for checking values from GT
!======================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkInterp : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(1,:) = ', TN_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(700,:) = ', TN_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(NPTS,:) = ', TN_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(1,:) = ', O_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(700,:) = ', O_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(NPTS,:) = ', O_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(1,:) = ', O2_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(700,:) = ', O2_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(NPTS,:) = ', O2_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(1,:) = ', N2_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(700,:) = ', N2_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(NPTS,:) = ', N2_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(1,:) = ', GLAt_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(700,:) = ', GLAt_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(NPTS,:) = ', GLAt_plasma_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(1,:) = ', GLOnd_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(700,:) = ', GLOnd_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(NPTS,:) = ', GLOnd_plasma_3d(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(1,:) = ', pz_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(700,:) = ', pz_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(NPTS,:) = ', pz_plasma_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(1,:) = ', V_east_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(700,:) = ', V_east_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(NPTS,:) = ', V_east_plasma(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(1,:) = ', V_south_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(700,:) = ', V_south_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(NPTS,:) = ', V_south_plasma(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(1,:) = ', V_upward_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(700,:) = ', V_upward_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(NPTS,:) = ', V_upward_plasma(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(1,:) = ', ilon1_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(700,:) = ', ilon1_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(NPTS,:) = ', ilon1_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(1,:) = ', ilon2_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(700,:) = ', ilon2_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(NPTS,:) = ', ilon2_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(1,:) = ', ilat1_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(700,:) = ', ilat1_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(NPTS,:) = ', ilat1_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(1,:) = ', ilat2_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(700,:) = ', ilat2_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(NPTS,:) = ', ilat2_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(1,:) = ', ispecial_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(700,:) = ', ispecial_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(NPTS,:) = ', ispecial_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(1,:) = ', ihl_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(700,:) = ', ihl_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(NPTS,:) = ', ihl_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(1,:) = ', ihu_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(700,:) = ', ihu_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(NPTS,:) = ', ihu_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, I3)' ) 'isFirstCallFixedHeight = ', isFirstCallFixedHeight


END SUBROUTINE checkFixedGeo


!====================================================================================================
END MODULE moduleDriverDebug
