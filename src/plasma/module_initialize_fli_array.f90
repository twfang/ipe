      module module_initialize_fli_array
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: initialize_fli_array
      contains
        subroutine initialize_fli_array ( )
          USE module_precision
          USE module_eldyn,ONLY: plas_fli !t,Je_3d
          USE module_physical_constants,ONLY:zero
          IMPLICIT NONE
!
          plas_fli(:,:,:,:) = zero
!t           IF ( sw_3DJ==1 )  Je_3d(IN:IS,mp,1:2) = zero
        end subroutine initialize_fli_array
      end module module_initialize_fli_array
