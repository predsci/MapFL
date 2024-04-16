c#######################################################################
      module magfld_func_def
c
c-----------------------------------------------------------------------
c ****** Definition of the analytic function that defines B.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Number of defined functions.
c
      integer, parameter :: number_of_functions=2
c
c ****** Mnemonics for defined functions.
c ****** There should be NUMBER_OF_FUNCTIONS of these.
c
      integer, parameter :: FUNC_TYPE_DIPOLE           =1
      integer, parameter :: FUNC_TYPE_PFSS_BKG         =2
c
      end module
c#######################################################################
      module magfld_func_index
c
c-----------------------------------------------------------------------
c ****** Index for the analytic function that defines B.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Selected function index.
c
      integer :: function_index=0
c
      end module
c#######################################################################
      module magfld_func_params
c
c-----------------------------------------------------------------------
c ****** Parameters for the analytic function that defines B.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Parameters for the DIPOLE function.
c
      real(r_typ) :: b0=1._r_typ
c
c ****** Parameters for the PFSS_BKG function.
c
      real(r_typ) :: mu
      real(r_typ) :: rss
c
      end module
c#######################################################################
      module magfld_func_namelist
c
c-----------------------------------------------------------------------
c ****** Namelist parameters for the analytic function that defines B.
c-----------------------------------------------------------------------
c
      use magfld_func_index
      use magfld_func_params
c
      implicit none
c
c ****** NAMELIST to read in the parameters.
c
      namelist /function_parameters/ function_index,
     &                               b0,
     &                               mu,rss
c
      end module
