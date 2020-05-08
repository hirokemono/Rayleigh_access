!>@file   rayleigh_restart_IO.f90
!!@brief  module rayleigh_restart_IO
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in 2004
!
!>@brief  Set file extension
!!
!!@verbatim
!!      character(len = kchara) function sel_rayleigh_file_name         &
!!     &                               (i_version, dir, int_id, postfix)
!!      subroutine sel_read_rayleigh_rst_params(dir, i_step, ra_rst)
!!        type(rayleigh_restart), intent(inout) :: ra_rst
!!          Version 0.99: "[dir]/[int_id]_[field_flag]"
!!          Version 1.x:  "[dir]/[int_id]/[field_flag]"
!!
!!      subroutine init_rayleigh_restart_input                          &
!!     &         (i_version, dir_name, i_step, fld_IO)
!!        type(field_IO), intent(inout) :: fld_IO
!!      subroutine set_rayleigh_rst_file_name(i_version, dir, i_step,   &
!!     &          field_name, iflag_ncomp, file_name)
!!@endverbatim
!
      module rayleigh_restart_IO
!
      use m_precision
!
      use t_field_data_IO
      use t_rayleigh_restart_IO
!      use t_base_field_labels
!      use t_explicit_term_labels
!
      implicit  none
!
!>      Number of components for scalar
      integer (kind = kint), parameter :: n_scalar = 1
!>      Number of components for vector
      integer (kind = kint), parameter :: n_vector = 3
!>      Number of components for solenoidal vector
      integer (kind = kint), parameter :: n_solenoid = 2
!
!>      Structure of field label
      type field_def
!>        Field name
        character(len=kchara) :: name
!>        Math expression
        character(len=kchara) :: math
!>        Number of component
        integer(kind = kint) ::  n_comp
      end type field_def
!
!>        Field label for velocity
!!         @f$ u_{i} @f$
      type(field_def), parameter :: velocity                            &
     &    = field_def(n_comp = n_vector,                                &
     &                name = 'velocity',                                &
     &                math = '$ u_{i} $')
!>        Field label for pressure @f$ p @f$
      type(field_def), parameter :: pressure                            &
     &    = field_def(n_comp = n_scalar,                                &
     &                name = 'pressure',                                &
     &                math = '$ u_{i} $')
!>        Field label for magnetic field
!!         @f$ B_{i} @f$
      type(field_def), parameter :: magnetic_field                      &
     &    = field_def(n_comp = n_vector,                                &
     &                name = 'magnetic_field',                          &
     &                math = '$ B_{i} $')
!>        Field label for temperature
!!         @f$ T @f$
      type(field_def), parameter :: temperature                         &
     &    = field_def(n_comp = n_scalar,                                &
     &                name = 'temperature',                             &
     &                math = '$ T $')
!
!>        Field label for explicit term for momentum at previous step
      type(field_def), parameter :: previous_momentum                   &
     &    = field_def(n_comp = n_vector,                                &
     &                name = 'previous_momentum',                       &
     &                math = '$ F_{u}^{N-1} $')
!>        Field label for explicit term for induction at previous step
      type(field_def), parameter :: previous_induction                  &
     &    = field_def(n_comp = n_vector,                                &
     &                name = 'previous_induction',                      &
     &                math = '$ F_{B}^{N-1} $')
!>        Field label for explicit term for heat at previous step
      type(field_def), parameter :: previous_heat                       &
     &    = field_def(n_comp = n_scalar,                                &
     &                name = 'previous_heat',                           &
     &                math = '$ F_{T}^{N-1} $')
!
      private :: set_rayleigh_restart_field
      private :: check_raylegh_rst_field
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      character(len = kchara) function sel_rayleigh_file_name           &
     &                               (i_version, dir, int_id, postfix)
!
      use rayleigh99_rst_param_IO
!
      integer(kind = kint), intent(in) :: i_version, int_id
      character(len=kchara), intent(in) :: dir, postfix
!
!
      if(i_version .lt. 1) then
        sel_rayleigh_file_name                                          &
     &      = set_rayleigh99_file_name(dir, int_id, postfix)
      else
        sel_rayleigh_file_name                                          &
     &      = set_rayleigh_file_name(dir, int_id, postfix)
      end if
!
      end function sel_rayleigh_file_name
!
!-----------------------------------------------------------------------
!
      subroutine sel_read_rayleigh_rst_params(dir, i_step, ra_rst)
!
      use rayleigh99_rst_param_IO
!
      integer(kind = kint), intent(in) :: i_step
      character(len = kchara), intent(in) :: dir
!
      type(rayleigh_restart), intent(inout) :: ra_rst
!
!
      if(ra_rst%i_version .lt. 1) then
        call read_rayleigh99_restart_params(dir, i_step, ra_rst)
      else
        call read_rayleigh_restart_params(dir, i_step, ra_rst)
      end if
!
      end subroutine sel_read_rayleigh_rst_params
!
!-----------------------------------------------------------------------
!
      subroutine init_rayleigh_restart_input                            &
     &         (i_version, dir_name, i_step, fld_IO)
!
      use calypso_mpi_4_rayleigh
!
      integer(kind = kint), intent(in) :: i_version, i_step
      character(len = kchara), intent(in) :: dir_name
!
      type(field_IO), intent(inout) :: fld_IO
!
      integer :: ilength
!
!
      if(my_rank .eq. 0) then
        call set_rayleigh_restart_field                                 &
     &     (i_version, dir_name, i_step, fld_IO)
      end if
!
      call MPI_Bcast(fld_IO%num_field_IO, 1,                            &
     &    CALYPSO_INTEGER, 0, CALYPSO_COMM, ierr_MPI)
      if(my_rank .ne. 0) call alloc_phys_name_IO(fld_IO)
!
      ilength = int(fld_IO%num_field_IO * kchara)
      call MPI_Bcast(fld_IO%fld_name, ilength,                          &
     &    CALYPSO_CHARACTER, 0, CALYPSO_COMM, ierr_MPI)
      call MPI_Bcast(fld_IO%num_comp_IO, fld_IO%num_field_IO,           &
     &    CALYPSO_INTEGER, 0, CALYPSO_COMM, ierr_MPI)
!
      call cal_istack_phys_comp_IO(fld_IO)
!      call alloc_phys_data_IO(fld_IO)
!
      end subroutine init_rayleigh_restart_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      logical function check_raylegh_rst_field                          &
     &                (i_version, dir, i_step, sclchar)
!
      integer(kind = kint), intent(in) :: i_version, i_step
      character(len = kchara), intent(in) :: dir
      character(len = kchara), intent(in) :: sclchar
!
      character(len = kchara) :: file1
!
!
      file1 = sel_rayleigh_file_name(i_version, dir, i_step, sclchar)
      check_raylegh_rst_field = check_file_exist(file1)
!
      end function check_raylegh_rst_field
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_rayleigh_restart_field                             &
     &         (i_version, dir, i_step, fld_IO)
!
      integer(kind = kint), intent(in) :: i_version, i_step
      character(len = kchara), intent(in) :: dir
!
      type(field_IO), intent(inout) :: fld_IO
!
!
      fld_IO%num_field_IO = 0
      call alloc_phys_name_IO(fld_IO)
!
      if(     check_raylegh_rst_field(i_version, dir, i_step, wchar)    &
     &  .and. check_raylegh_rst_field(i_version, dir, i_step, zchar))   &
     &   call append_phys_name_IO(velocity%name, n_solenoid, fld_IO)
      if(     check_raylegh_rst_field(i_version, dir, i_step, cchar)    &
     &  .and. check_raylegh_rst_field(i_version, dir, i_step, achar))   &
     &   call append_phys_name_IO(magnetic_field%name, n_solenoid,      &
     &                            fld_IO)
      if(check_raylegh_rst_field(i_version, dir, i_step, pchar))        &
     &   call append_phys_name_IO(pressure%name, n_scalar, fld_IO)
      if(check_raylegh_rst_field(i_version, dir, i_step, tchar))        &
     &   call append_phys_name_IO(temperature%name, n_scalar, fld_IO)
!
!
      if(     check_raylegh_rst_field(i_version, dir, i_step, wabchar)  &
     &  .and. check_raylegh_rst_field(i_version, dir, i_step, zabchar)) &
     &    call append_phys_name_IO(previous_momentum%name, n_solenoid,  &
     &                             fld_IO)
      if(     check_raylegh_rst_field(i_version, dir, i_step, cabchar)  &
     &  .and. check_raylegh_rst_field(i_version, dir, i_step, aabchar)) &
     &    call append_phys_name_IO(previous_induction%name, n_solenoid, &
     &                             fld_IO)
!      if(check_raylegh_rst_field(i_version, dir, i_step, pabchar))     &
!     &   call append_phys_name_IO(previous_pressure%name, n_scalar,    &
!     &                            fld_IO)
      if(check_raylegh_rst_field(i_version, dir, i_step, tabchar))      &
     &   call append_phys_name_IO(previous_heat%name, n_scalar, fld_IO)
!
      end subroutine set_rayleigh_restart_field
!
!-----------------------------------------------------------------------
!
      subroutine set_rayleigh_rst_file_name(i_version, dir, i_step,     &
     &          field_name, iflag_ncomp, file_name)
!
      integer(kind = kint), intent(in) :: i_version, i_step
      character(len = kchara), intent(in) :: dir
      character(len = kchara), intent(in) :: field_name
!
      integer(kind = kint), intent(inout) :: iflag_ncomp
      character(len = kchara), intent(inout) :: file_name(2)
!
      if(field_name .eq. velocity%name) then
        iflag_ncomp = 2
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         wchar)
        file_name(2) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         zchar)
      else if(field_name .eq. pressure%name) then
        iflag_ncomp = 1
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         pchar)
      else if(field_name .eq. temperature%name) then
        iflag_ncomp = 1
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         tchar)
      else if(field_name .eq. magnetic_field%name) then
        iflag_ncomp = 2
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         cchar)
        file_name(2) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         achar)
!
      else if(field_name .eq. previous_momentum%name) then
        iflag_ncomp = 2
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         wabchar)
        file_name(2) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         zabchar)
!      else if(field_name .eq. previous_pressure%name) then
!        iflag_ncomp = 1
!        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step, &
!     &                                         cchar)
      else if(field_name .eq. previous_heat%name) then
        iflag_ncomp = 1
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         tabchar)
      else if(field_name .eq. previous_induction%name) then
        iflag_ncomp = 2
        file_name(1) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         cabchar)
        file_name(2) =  sel_rayleigh_file_name(i_version, dir, i_step,  &
     &                                         aabchar)
      end if
!
      end subroutine set_rayleigh_rst_file_name
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      logical function check_file_exist(file_name)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), parameter :: id_file = 255
!
!
      open (id_file, file=file_name, status='old', err=99)
      close(id_file)
!
      check_file_exist = .TRUE.
      return
!
  99  continue
      check_file_exist = .FALSE.
      return
!
      end function check_file_exist
!
!  ---------------------------------------------------------------------
!
      end module rayleigh_restart_IO
