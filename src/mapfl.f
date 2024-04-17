c#######################################################################
c
c-----------------------------------------------------------------------
c
c     __  __      _      ____    _____   _
c    |  \/  |    / \    |  _ \  |  ___| | |
c    | |\/| |   / _ \   | |_) | | |_    | |
c    | |  | |  / ___ \  |  __/  |  _|   | |___
c    |_|  |_| /_/   \_\ |_|     |_|     |_____|
c
c
c ****** MAPFL: Map field lines through a 3D vector field.
c
c     Authors:  Zoran Mikic
c               Jon Linker
c               Cooper Downs
c               Roberto Lionello
c               Ronald M. Caplan
c               Emily Mason
c
c     Predictive Science Inc.
c     www.predsci.com
c     San Diego, California, USA 92121
c
c#######################################################################
c Copyright 2024 Predictive Science Inc.
c
c Licensed under the Apache License, Version 2.0 (the "License");
c you may not use this file except in compliance with the License.
c You may obtain a copy of the License at
c
c    http://www.apache.org/licenses/LICENSE-2.0
c
c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
c implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c#######################################################################
c
c#######################################################################
      module ident
c
c-----------------------------------------------------------------------
c ****** Set the name, version, and date of code.
c-----------------------------------------------------------------------
c
c
      character(*), parameter :: cname='MAPFL'
      character(*), parameter :: cvers='2.1.0'
      character(*), parameter :: cdate='04/16/2024'
c
      end module
c#######################################################################
      module number_types
c
c-----------------------------------------------------------------------
c ****** Set precisions for REALs.
c-----------------------------------------------------------------------
c
      use iso_fortran_env
c
c ****** Use double precision.
c
      integer, parameter :: r_typ = REAL64
c
      end module
c#######################################################################
      module spline_def
c
c-----------------------------------------------------------------------
c ****** Definition of cubic spline data structures.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ***** 1D spline structure.
c
      type :: spl1d
        integer :: nx
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: f
        real(r_typ), dimension(:), pointer :: fxx
      end type
c
c ***** 2D spline structure.
c
      type :: spl2d
        integer :: nx
        integer :: ny
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:,:), pointer :: f
        real(r_typ), dimension(:,:), pointer :: fxx
        real(r_typ), dimension(:,:), pointer :: fyy
        real(r_typ), dimension(:,:), pointer :: fxxyy
      end type
c
c ***** 3D spline structure.
c
      type :: spl3d
        integer :: nx
        integer :: ny
        integer :: nz
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:), pointer :: z
        real(r_typ), dimension(:,:,:), pointer :: f
        real(r_typ), dimension(:,:,:), pointer :: fxx
        real(r_typ), dimension(:,:,:), pointer :: fyy
        real(r_typ), dimension(:,:,:), pointer :: fzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyy
        real(r_typ), dimension(:,:,:), pointer :: fxxzz
        real(r_typ), dimension(:,:,:), pointer :: fyyzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyyzz
      end type
c
      end module
c#######################################################################
      module invint_def
c
c-----------------------------------------------------------------------
c ****** Definition of an inverse interpolation table data structure.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      type :: itab
        integer :: n
        real(r_typ), dimension(:), pointer :: f
        real(r_typ) :: d
      end type
c
      end module
c#######################################################################
      module locate_interval_interface
      interface
        function locate_interval (n,x,xv,tab,ierr)
        use number_types
        use invint_def
        implicit none
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        type(itab), optional :: tab
        integer, optional :: ierr
        integer :: locate_interval
        intent(in) :: n,x,xv,tab
        end
      end interface
      end module
c#######################################################################
      module evaluate_spline_3d_interface
      interface
        function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
        use number_types
        use spline_def
        use invint_def
        use locate_interval_interface
        type(spl3d) :: s
        real(r_typ) :: x,y,z
        type(itab), optional :: tabx,taby,tabz
        real(r_typ) :: evaluate_spline_3d
        end
      end interface
      end module
c#######################################################################
      module debug
c
      implicit none
c
c ****** Debugging level.
c
      integer :: debug_level=0
c
      end module
c#######################################################################
      module constants
c
c-----------------------------------------------------------------------
c ****** Constants.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: halfpi=.5_r_typ*pi
      real(r_typ), parameter :: twopi=2._r_typ*pi
c
      end module
c#######################################################################
      module types
c
c-----------------------------------------------------------------------
c ****** Definition of data structures.
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use invint_def
      use spline_def
c
      implicit none
c
c ****** Maximum number of dimensions.
c
      integer, parameter, private :: ndim_max=3
c
c ****** Inverse interpolation table structure for a vector field.
c
      type :: vtab
        type(itab), dimension(ndim_max) :: c
      end type
c
c ****** Vector spline structure.
c
      type :: vspl3d
        type(spl3d) :: r
        type(spl3d) :: t
        type(spl3d) :: p
      end type
c
c ****** Magnetic field vector structure.
c
      type :: vec
        type(sds) :: r
        type(sds) :: t
        type(sds) :: p
        integer :: nrs
        integer :: nts
        integer :: nps
        real(r_typ), dimension(:), pointer :: rs
        real(r_typ), dimension(:), pointer :: ts
        real(r_typ), dimension(:), pointer :: ps
        real(r_typ) :: lim0(ndim_max)
        real(r_typ) :: lim1(ndim_max)
        type(vtab), dimension(ndim_max) :: inv
        real(r_typ), dimension(:), pointer :: drs
        real(r_typ), dimension(:), pointer :: dts
        real(r_typ), dimension(:), pointer :: dps
        real(r_typ), dimension(:), pointer :: sts
        type(itab) :: rs_invtab
        type(itab) :: ts_invtab
        type(itab) :: ps_invtab
        logical :: cubic=.false.
        type(vspl3d) :: spl
        logical :: b_is_32bit
      end type
c
c ****** Data structure to hold file names of a vector component.
c
      type :: vfile
        character(512) :: r=' '
        character(512) :: t=' '
        character(512) :: p=' '
      end type
c
c ****** Initial size for the field line trace buffer.
c
      integer, parameter, private :: fl_buffer_size=1000
c
c ****** Trajectory structure definition.
c
      type :: traj
        integer :: ndim
        integer :: initial_size=fl_buffer_size
        integer :: size
        integer :: npts
        type(rp1d), dimension(:), pointer :: x
      end type
c
c ****** Dual representation Cartesian and spherical position vector.
c
      type :: csvec
        real(r_typ), dimension(3) :: c
        real(r_typ), dimension(3) :: s
      end type
c
c ****** "Inside domain" structure.
c
      type :: inout
        logical :: domain
        logical :: r0
        logical :: r1
        logical :: t0
        logical :: t1
        logical :: p0
        logical :: p1
        logical :: r
        logical :: t
        logical :: p
      end type
c
c ****** Field line integration parameters.
c
      type :: flparam
        logical :: variable=.true.
        integer :: direction
        logical :: direction_is_along_b
        real(r_typ) :: min=0.001_r_typ
        real(r_typ) :: max=0.1_r_typ
        real(r_typ) :: over_rc=0.0025_r_typ
        real(r_typ) :: lmax=100.0_r_typ
        logical :: limit_by_local_mesh=.true.
        real(r_typ) :: local_mesh_factor=1.0_r_typ
        real(r_typ) :: max_increase_factor
        real(r_typ) :: max_decrease_factor
        integer :: short_fl_min_points
        integer :: short_fl_max_tries
        real(r_typ) :: short_fl_shrink_factor
        real(r_typ) :: predictor_min_clip_fraction
      end type
c
c ****** Preferences for Python/f2py (used in other source files).
c
      type :: preferences
        logical :: cubic_vec_field
        logical :: var_dstep
        real(r_typ) :: dstep
        logical :: auto_minmax_dstep
        real(r_typ) :: min_dstep
        real(r_typ) :: max_dstep
        real(r_typ) :: dstep_mult
        logical :: limit_by_local_mesh
        real(r_typ) :: local_mesh_factor
        real(r_typ) :: max_length
        logical :: direction_along_vec_field
        logical :: trace_from_slice_forward
        logical :: trace_from_slice_backward
      end type
c
      end module
c#######################################################################
      module mesh
c
c-----------------------------------------------------------------------
c ****** Meshes.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Secondary (r,t,p) meshes.
c
      integer :: nrss=54
      integer :: ntss=181
      integer :: npss=361
c
      real(r_typ), dimension(:), pointer :: rss
      real(r_typ), dimension(:), pointer :: tss
      real(r_typ), dimension(:), pointer :: pss
c
      end module
c#######################################################################
      module field
c
c-----------------------------------------------------------------------
c ****** Magnetic field storage.
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
c ****** Structure that holds the magnetic field.
c
      type(vec) :: b
c
      end module
c#######################################################################
      module vars
c
c-----------------------------------------------------------------------
c ****** Input variables, switches, etc.
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use types
c
      implicit none
c
c ****** Spherical geometry domain limits.
c
      real(r_typ) :: domain_r_min=1._r_typ
      real(r_typ) :: domain_r_max=300._r_typ
c
c ****** New output mesh flags.
c
      logical :: new_r_mesh=.false.
      logical :: new_t_mesh=.false.
      logical :: new_p_mesh=.false.
c
c ****** New output mesh limits.
c
      real(r_typ) :: r0=0.
      real(r_typ) :: r1=0.
      real(r_typ) :: t0=0.
      real(r_typ) :: t1=0.
      real(r_typ) :: p0=0.
      real(r_typ) :: p1=0.
c
c ****** Flag to use tri-cubic interpolation (when .TRUE.) or
c ****** simple linear interpolation (when .FALSE.) to
c ****** interpolate B between mesh points.
c
      logical :: cubic=.false.
c
c ****** Field line integration.
c
      type(flparam) :: ds
      logical :: set_ds_automatically=.true.
      real(r_typ) :: dsmult=1.0_r_typ
c
c ****** Flag to request a mapping on a slice.
c
      logical :: trace_slice=.false.
c
c ****** Parameters for the slice mapping.
c
      logical :: slice_coords_are_xyz=.false.
c
c ****** Flag to specify whether the tracing direction is along B
c ****** or along the direction of increasing radius.
c
      logical :: trace_slice_direction_is_along_b=.true.
c
c ****** Flags to request tracing in the forward and backward
c ****** directions.
c
      logical :: trace_from_slice_forward=.false.
      logical :: trace_from_slice_backward=.false.
c
c ****** Names of the slice coodrinates.
c
      character, dimension(3) :: slice_coord_name
c
c ****** Structures that hold the slice coordinates.
c
      type (sds) :: slice_c1,slice_c2,slice_c3
c
c ****** Increment to compute Q directly on a slice.
c
      real(r_typ) :: q_increment_h=0.0001_r_typ
c
c ****** Flag to request a 3D mapping.
c
      logical :: trace_3d=.false.
c
c ****** Switch to use an analytic function to define the magnetic
c ****** field.
c
      logical :: use_analytic_function=.false.
c
c ****** Flag to use 32-bit HDF output files.
c
      logical :: hdf32=.true.
c
c ****** Flag to write field line traces originating from a slice
c ****** to individual HDF output files.
c
      logical :: write_traces_to_hdf=.false.
c
c ****** String used for the root file name for the field line
c ****** traces.
c
      character(64) :: write_traces_root='fl'
c
c ****** Flag to write the Cartesian (x,y,z) coordinates for
c ****** field line traces.
c
      logical :: write_traces_as_xyz=.false.
c
c ***** File type for field line tracing output.
c
      character(3) :: fmt='hdf'
c
c ***** Number of segnments along which to search for dips.
c
      integer :: ns_dips=10
c
c ****** Integrate scalar field along field line.
c
      logical :: integrate_along_fl=.false.
c
      end module
c#######################################################################
      module diags
c
c-----------------------------------------------------------------------
c ****** Variables that control diagnostic output.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Number of iterations between prints of diagnostics
c ****** during execution.
c
      integer :: diagnostic_interval=1000
c
      end module
c#######################################################################
      module files
c
c-----------------------------------------------------------------------
c ****** File names.
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
      type(vfile) :: bfile
c
      character(512) :: rffile=' ',tffile=' ',pffile=' ',effile=' '
      character(512) :: kffile=' ',qffile=' ',slogqffile=' '
      character(512) :: rbfile=' ',tbfile=' ',pbfile=' ',ebfile=' '
      character(512) :: kbfile=' ',qbfile=' ',slogqbfile=' '
      character(512) :: lffile=' ',lbfile=' '
c
c ****** File names for the r, t, and p meshes.
c
      character(512) :: mesh_file_r=' '
      character(512) :: mesh_file_t=' '
      character(512) :: mesh_file_p=' '
c
      type(vfile) :: volume3d_output_file
c
      type(vfile) :: slice_input_file
      type(vfile) :: slice_output_file_forward
      type(vfile) :: slice_output_file_backward
c
c ****** File name for the analytic function parameters.
c
      character(512) :: function_params_file
     &                                  ='magfield_function_params.dat'
c
c ****** File names for slice output quantities.
c
      character(512) :: slice_q_output_file=' '
      character(512) :: slice_length_output_file=' '
c
c ****** File name for the output coronal hole map.
c
      character(512) :: ch_map_output_file='ch.h5'
c
c ****** File name for the output 3D coronal hole map.
c
      character(512) :: ch_map_3d_output_file='ch3d.h5'
c
c ****** File name for the output 3D dips map.
c
      character(512) :: dips_map_3d_output_file='dips3d.h5'
c
c ****** File name for the scalar field to be integrated along fl.
c
      character(512) :: scalar_input_file=' '
c
      end module
c#######################################################################
      module field_line_params
c
c-----------------------------------------------------------------------
c ****** Parameters that control the field line integration.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Parameters that control variable step-size tracing.
c-----------------------------------------------------------------------
c
c ****** These factors control how much the field line integration
c ****** step size can change from one step to another for the
c ****** case when a variable step size is being used.
c
c ****** MAX_INCREASE_FACTOR should be greater than 1, and
c ****** MAX_DECREASE_FACTOR should be less than 1.

      real(r_typ), parameter :: max_increase_factor=1.5_r_typ
      real(r_typ), parameter :: max_decrease_factor=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Parameters that control tracing of short field lines.
c-----------------------------------------------------------------------
c
c ****** If a field line trace has a smaller number of points
c ****** than SHORT_FL_MIN_POINTS, the integration step size
c ****** is decreased by the factor SHORT_FL_SHRINK_FACTOR,
c ****** and it is retraced, up to a maximum number of tries
c ****** equal to SHORT_FL_MAX_TRIES.
c
c ****** SHORT_FL_SHRINK_FACTOR should be less than 1.
c
      integer, parameter :: short_fl_min_points=10
      integer, parameter :: short_fl_max_tries=5
      real(r_typ), parameter :: short_fl_shrink_factor=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Parameters that control clipping to boundaries.
c-----------------------------------------------------------------------
c
c ****** The factor PREDICTOR_MIN_CLIP_FRACTION determines when to
c ****** clip a trace to the radial boundary in the predictor.
c ****** When the normalized distance to the r boundary (as a
c ****** fraction of the current step size) is less than
c ****** PREDICTOR_MIN_CLIP_FRACTION, the field line is clipped
c ****** to the boundary in the predictor without doing a
c ****** corrector step.  This number should be between 0 and 1.
c
      real(r_typ), parameter :: predictor_min_clip_fraction=.1_r_typ
c
c-----------------------------------------------------------------------
c ****** Maximum number of "bad" field line traces after
c ****** which to terminate.  Set to -1 to disable the termination.
c-----------------------------------------------------------------------
c
      integer, parameter :: max_bad_fieldlines=-1
c
      end module
c#######################################################################
      module step_size_stats
c
c-----------------------------------------------------------------------
c ****** Variable step size statistics.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      integer :: gather_stats = 0
c
      integer(8) :: stat_n=0
      real(r_typ) :: stat_ds_sum=0._r_typ
      real(r_typ) :: stat_ds_avg=0._r_typ
      real(r_typ) :: stat_ds_min=huge(0._r_typ)
      real(r_typ) :: stat_ds_max=0._r_typ
c
      end module
c#######################################################################
      module openmp_vars
c
c-----------------------------------------------------------------------
c ****** Variables to control OpenMP parallelization.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Number of iterations to do in each thread.
c
      integer :: iterations_per_thread=500
c
      end module
c#######################################################################
      module params
c
c-----------------------------------------------------------------------
c ****** Parameters.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      character(len=:), allocatable :: infile
      integer :: verbose = 0
c
      end module
c#######################################################################
      module interp_interface
      interface
        subroutine interp (n,x,xv,i,ip1,alpha,tab)
        use number_types
        use invint_def
        use locate_interval_interface
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        integer :: i
        integer :: ip1
        real(r_typ) :: alpha
        type(itab), optional :: tab
        intent(in) :: n,x,xv,tab
        intent(out) :: i,ip1,alpha
        end
      end interface
      end module
c#######################################################################
      module tracefl_interface
      interface
        subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,
     &                      traced_to_r_boundary,xt)
        use number_types
        use types
        type(vec) :: b
        type(flparam) :: ds
        real(r_typ), dimension(3) :: s0,s1
        real(r_typ), dimension(3) :: bs0,bs1
        real(r_typ) :: s
        logical :: traced_to_r_boundary
        type(traj), optional :: xt
        intent(in) :: b,ds,s0
        intent(out) :: s1,bs0,bs1,s,traced_to_r_boundary
        end
      end interface
      end module
c#######################################################################
      module integrate_fl
c
c-----------------------------------------------------------------------
c ****** Internal variables for integration along field line
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
      logical :: do_integral_along_fl=.false.
      type(sds) :: scalar_field
      type(vtab) :: inv_sf

c
      end module
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
      program MAPFL
c
c-----------------------------------------------------------------------
c
      use ident
      use params
      use types
      use files
      use mesh
      use field
      use vars
      use field_line_params
      use step_size_stats
      use debug
      use magfld_func_def
      use magfld_func_index
      use magfld_func_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i
      real(r_typ) :: ch_map_r=1.0_r_typ
      logical :: trace_fwd=.false.
      logical :: trace_bwd=.false.
      logical :: compute_q_on_slice=.false.
      logical :: compute_ch_map=.false.
      logical :: compute_ch_map_3d=.false.
      logical :: compute_dips_map_3d=.false.
      character(256) :: errline=' '
      integer arglen
c
c-----------------------------------------------------------------------
c
      namelist /datum/
     &  debug_level, use_analytic_function, function_params_file,
     &  domain_r_min, domain_r_max, bfile,
     &  cubic, ds, set_ds_automatically,
     &  dsmult, trace_fwd, trace_bwd,
     &  rffile, tffile, pffile, effile, kffile, qffile, lffile,
     &  rbfile, tbfile, pbfile, ebfile, kbfile, qbfile, lbfile,
     &  new_r_mesh, mesh_file_r, nrss, r0,r1, new_t_mesh, mesh_file_t,
     &  ntss, t0,t1, new_p_mesh, mesh_file_p, npss, p0,p1, trace_3d,
     &  volume3d_output_file, trace_slice, slice_coords_are_xyz,
     &  trace_slice_direction_is_along_b, compute_q_on_slice,
     &  q_increment_h, slice_input_file, trace_from_slice_forward,
     &  slice_output_file_forward, trace_from_slice_backward,
     &  slice_output_file_backward, slice_q_output_file,
     &  slice_length_output_file, compute_ch_map, ch_map_r,
     &  ch_map_output_file, compute_ch_map_3d, ch_map_3d_output_file,
     &  write_traces_to_hdf, write_traces_root, write_traces_as_xyz,
     &  compute_dips_map_3d, dips_map_3d_output_file, ns_dips,
     &  slogqffile,slogqbfile,integrate_along_fl,scalar_input_file,
     &  verbose, function_index, b0, mu
c
c-----------------------------------------------------------------------
c
c ****** Get input filename from command line.
c
      if (command_argument_count().ge.1) then
        call GET_COMMAND_ARGUMENT(1,length=arglen)
        allocate(character(arglen) :: infile)
        call GET_COMMAND_ARGUMENT(1,value=infile)
      else
        allocate(character(9):: infile)
        infile='mapfl.in'
      end if
c
      call ffopen (8,infile,'r',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in MAPFL:'
        write (*,*) '### The input file does not exist'//
     &              ' or cannot be read.'
        write (*,*) 'File name: ',trim(infile)
        call exit (1)
      end if
c
c ****** Read the input file.
c
      call ffopen (8,trim(infile),'r',ierr)
      read(8,datum,iostat=ierr)
      if (ierr.ne.0) then
        backspace (8)
        read (1,fmt='(A)') errline
        write (*,*)
        write (*,*) '### ERROR reading input file:'
        write (*,*) '### The following line has a problem:'
        write (*,*)
        write (*,*) trim(errline)
        write (*,*)
        write (*,*) '###'
        call exit (1)
      endif
      write (*,*)
      write (*,*) '### Input file contents:'
      write (*,*)
      write(*,datum)
      close (8)
c
c
c ****** Write the NAMELIST parameter values to file.
c
      call ffopen (8,'mapfl_parameters_used.out','rw',ierr)
      write (8,datum)
      close (8)
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
      end if
c
c ****** Read the parameters that define the analytic magnetic
c ****** field function, if requested.
c
      if (use_analytic_function) then
        call process_function_params
      end if
c
c ****** Set the field line integration parameters.
c
      ds%max_increase_factor=max_increase_factor
      ds%max_decrease_factor=max_decrease_factor
      ds%predictor_min_clip_fraction=predictor_min_clip_fraction
      ds%short_fl_min_points=short_fl_min_points
      ds%short_fl_max_tries=short_fl_max_tries
      ds%short_fl_shrink_factor=short_fl_shrink_factor
c
c ****** Read the magnetic field.
c
      if (.not.use_analytic_function) then
        call readb (bfile,b)
      end if
c
c ****** Set the trace output format based on input br
c ****** (for analytic function, sets to hdf)
c
      i=index(bfile%r,'.h');
      if (bfile%r(i+1:i+2).eq.'h5') then
        fmt='h5'
      endif
c
c ****** Set the radial domain limits to those specified.
c
      b%lim0(1)=max(b%lim0(1),domain_r_min)
      b%lim1(1)=min(b%lim1(1),domain_r_max)
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Domain limits:'
        write (*,*) 'Lower boundary value: ',b%lim0(1)
        write (*,*) 'Upper boundary value: ',b%lim1(1)
      end if
c
c ****** Make the new r, t, and p meshes.
c
      call make_new_meshes (b)
c
c ****** Set the default step size.
c
      call set_ds (b)
c
c ****** Set the flag to gather step size statistics.
c
      gather_stats=verbose
c
c ****** Setup the field to integrate along if requested.
c
      if (integrate_along_fl) call set_up_integration
c
c ****** Trace the field lines forward, if requested.
c

      if (trace_fwd) call map_forward
c
c ****** Trace the field lines backward, if requested.
c
      if (trace_bwd) call map_backward
c
c ****** Map the field lines from a 3D rectilinear volume,
c ****** if requested.
c
      if (trace_3d) call map_3d
c
c ****** Map the field lines from a slice, if requested,
c ****** or determine Q on the slice, if requested.
c
      if (trace_slice) then
        call read_slice_coordinates
        if (compute_q_on_slice) then
          call get_q_on_slice
        end if
        call map_slice
        call deallocate_slice_coordinates
      end if
c
c ****** Compute a coronal hole map, if requested.
c
      if (compute_ch_map) then
        call get_ch_map (ch_map_r)
      end if
c
c ****** Compute a 3D coronal hole map, if requested.
c
      if (compute_ch_map_3d) then
        call get_ch_map_3d
      end if
c
c ****** Compute a 3D dips map, if requested.
c
      if (compute_dips_map_3d) then
        call get_dips_map_3d
      end if
c
      if (verbose.gt.0) then
        stat_ds_avg=0.
        if (stat_n.ne.0) stat_ds_avg=stat_ds_sum/stat_n
        write (*,*)
        write (*,*) '### Field line integration step size statistics:'
        write (*,*) 'Number of field line segments = ',stat_n
        write (*,*) 'Minimum step size used = ',stat_ds_min
        write (*,*) 'Maximum step size used = ',stat_ds_max
        write (*,*) 'Average step size used = ',stat_ds_avg
      end if
c
      call exit (0)
c
      end
c#######################################################################
      subroutine process_function_params
c
c-----------------------------------------------------------------------
c
c ****** Process parameters for the analytic magnetic field
c ****** function.
c
c-----------------------------------------------------------------------
c
      use number_types
      use ident
      use params
      use vars
      use files
      use field
      use magfld_func_def
      use magfld_func_index
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Using an analytic magnetic field'//
     &              ' function.'
      end if
c
c ****** Check that the requested function index is valid.
c
      if (function_index.lt.1.or.
     &    function_index.gt.number_of_functions) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### An invalid function index was requested:'
        write (*,*) 'FUNCTION_INDEX = ',function_index
        call exit (1)
      end if
c
      if (verbose.gt.0) then
        write (*,*) '### Using the analytic magnetic field'//
     &              ' function with index = ',function_index
      end if
c
c ****** Check that the selected geometry is consistent with the
c ****** use of an analytic field.
c
      if (.not.(new_r_mesh.and.new_t_mesh.and.new_p_mesh)) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must specify the r, t, and p meshes to use.'
        call exit (1)
      end if
c
      if (ds%limit_by_local_mesh) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must not attempt to limit the integration'//
     &              ' step size by the local'
        write (*,*) 'magnetic field mesh.'
        call exit (1)
      end if
c
      if (set_ds_automatically) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must not attempt to set the integration'//
     &              ' step size automatically.'
        call exit (1)
      end if
c
c ****** Load the domain limits into the B structure.
c
      b%lim0(1)=domain_r_min
      b%lim1(1)=domain_r_max
      b%lim0(2)=0.
      b%lim1(2)=pi
      b%lim0(3)=0.
      b%lim1(3)=twopi
c
      return
      end
c#######################################################################
      subroutine readb (bfile,b)
c
c-----------------------------------------------------------------------
c
c ****** Read the magnetic field from the files specified by
c ****** BFILE into the magnetic field vector structure B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use vars
      use params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vfile) :: bfile
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Read the magnetic field components.
c
c ****** Br.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Reading data file: ',trim(bfile%r)
      end if
c
      call rdhdf (bfile%r,b%r,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Br.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
c
      if (b%r%ndim.ne.3.or..not.b%r%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Br file.'
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
c
c ****** Bt.
c
      if (verbose.gt.0) then
        write (*,*) 'Reading data file: ',trim(bfile%t)
      end if
c
      call rdhdf (bfile%t,b%t,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bt.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
c
      if (b%t%ndim.ne.3.or..not.b%t%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bt file.'
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
c
c ****** Bp.
c
      if (verbose.gt.0) then
        write (*,*) 'Reading data file: ',trim(bfile%p)
      end if
c
      call rdhdf (bfile%p,b%p,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bp.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
c
      if (b%p%ndim.ne.3.or..not.b%p%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bp file.'
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
c
c ****** Set the type of magnetic field files.
      if (verbose.gt.0) then
        write (*,*) 'Setting btype.'
      end if
c
      call set_btype (b)
c
c ****** Build the inverse interpolation tables.
c

      if (verbose.gt.0) then
        write (*,*) 'Building inverse tables.'
      end if
c
      call build_inverse_tables (b%r,b%inv(1))
      call build_inverse_tables (b%t,b%inv(2))
      call build_inverse_tables (b%p,b%inv(3))
c
c ****** If cubic spline interpolation was requested, get the
c ****** spline coefficients.
c
      if (cubic) then
        b%cubic=.true.
        if (verbose.gt.0) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Br ...'
        end if
        call compute_spline_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),
     &                          b%r%scales(1)%f,
     &                          b%r%scales(2)%f,
     &                          b%r%scales(3)%f,
     &                          b%r%f,b%spl%r)
        if (verbose.gt.0) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Bt ...'
        end if
        call compute_spline_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),
     &                          b%t%scales(1)%f,
     &                          b%t%scales(2)%f,
     &                          b%t%scales(3)%f,
     &                          b%t%f,b%spl%t)
        if (verbose.gt.0) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//
     &                ' for Bp ...'
        end if
        call compute_spline_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),
     &                          b%p%scales(1)%f,
     &                          b%p%scales(2)%f,
     &                          b%p%scales(3)%f,
     &                          b%p%f,b%spl%p)
      else
        b%cubic=.false.
      end if
c
      return
      end
c#######################################################################
      subroutine set_btype (b)
c
c-----------------------------------------------------------------------
c
c ****** Determine the primary (r,t,p) scales and the mesh limits
c ****** from the type of magnetic field in structure B, and
c ****** store them in structure B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use constants
      use params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
c ****** Tolerance for checking the bounds of t and p scales.
c ****** This should be set to several times the roundoff in
c ****** 32-bit representations of pi.
c
      real(r_typ), parameter :: eps=2.e-6_r_typ
c
c ****** The values of pi and 2*pi using 32-bit precision.
c
      real(REAL32), parameter :: pi_r4=pi
      real(REAL32), parameter :: twopi_r4=twopi
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: z
      integer :: n1,n2,n3
      logical :: add_phi_point
c
c-----------------------------------------------------------------------
c
c ****** Check the type of magnetic field files read in.
c
      if (verbose.gt.0) then
        write (*,*) 'SET_BTYPE: Checking type of magnetic field files.'
      end if
c
c ****** Check for new MAS code files.
c
      if (b%r%dims(1).eq.b%t%dims(1)+1.and.
     &    b%r%dims(1).eq.b%p%dims(1)+1.and.
     &    b%r%dims(2).eq.b%t%dims(2)-1.and.
     &    b%r%dims(2).eq.b%p%dims(2)  .and.
     &    b%r%dims(3).eq.b%t%dims(3)  .and.
     &    b%r%dims(3).eq.b%p%dims(3)-1) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected new MAS code type.'
        end if
c
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for old MAS code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)+1.and.
     &         b%r%dims(1).eq.b%p%dims(1)+1.and.
     &         b%r%dims(2).eq.b%t%dims(2)-1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected old MAS code type.'
        end if
c
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
c ****** Do not add a phi point if the phi interval already includes
c ****** the whole interval. Can occur if grid modified outside mapfl.
c
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose.gt.0) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
c
c ****** Check for new POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)+1) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected POT3D code type.'
        end if
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)-1
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%p%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for new-old POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)-1) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected new-old POT3D type.'
        end if
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
        add_phi_point=.false.
c
c ****** Check for old POT3D code files.
c
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.
     &         b%r%dims(2).eq.b%p%dims(2)  .and.
     &         b%r%dims(3).eq.b%t%dims(3)  .and.
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected old POT3D code type.'
        end if
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
c ****** Do not add a phi point if the phi interval already includes
c ****** the whole interval. Can occur if grid modified outside mapfl.
c
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose.gt.0) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
c
c ****** Check for non-staggered files.
c
      else if (b%r%dims(1).eq.b%t%dims(1).and.
     &         b%r%dims(1).eq.b%p%dims(1).and.
     &         b%r%dims(2).eq.b%t%dims(2).and.
     &         b%r%dims(2).eq.b%p%dims(2).and.
     &         b%r%dims(3).eq.b%t%dims(3).and.
     &         b%r%dims(3).eq.b%p%dims(3)) then
c
        if (verbose.gt.0) then
          write(*,*) 'SET_BTYPE: Detected non-staggered type.'
        end if
c
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
c
        b%rs=>b%r%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
c
c ****** Do not add a phi point if the phi interval already includes
c ****** the whole interval. Can occur if grid modified outside mapfl.
c
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose.gt.0) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
c
      else
c
c ****** Invalid file type.
c
        write (*,*)
        write (*,*) '### ERROR in SET_BTYPE:'
        write (*,*) '### Unrecognized magnetic field file staggering.'
        write (*,*) '  br resolution:', b%r%dims
        write (*,*) '  bt resolution:', b%t%dims
        write (*,*) '  bp resolution:', b%p%dims
        call exit (1)
c
      end if
c
c ****** If appropriate, add a point in the phi dimension to
c ****** take care of periodic wrap-around.
c
      if (add_phi_point) then
c
        if (verbose.gt.0) then
          write (*,*) 'SET_BTYPE: Adding phi point.'
        end if
c
        n1=b%r%dims(1)
        n2=b%r%dims(2)
        n3=b%r%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%r%f(:,:,:)
        f(:,:,n3+1)=b%r%f(:,:,1)
        z(1:n3)=b%r%scales(3)%f(:)
        z(n3+1)=b%r%scales(3)%f(1)+twopi
        deallocate (b%r%f)
        deallocate (b%r%scales(3)%f)
        b%r%dims(3)=n3+1
        b%r%f=>f
        b%r%scales(3)%f=>z
c
        n1=b%t%dims(1)
        n2=b%t%dims(2)
        n3=b%t%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%t%f(:,:,:)
        f(:,:,n3+1)=b%t%f(:,:,1)
        z(1:n3)=b%t%scales(3)%f(:)
        z(n3+1)=b%t%scales(3)%f(1)+twopi
        deallocate (b%t%f)
        deallocate (b%t%scales(3)%f)
        b%t%dims(3)=n3+1
        b%t%f=>f
        b%t%scales(3)%f=>z
c
        n1=b%p%dims(1)
        n2=b%p%dims(2)
        n3=b%p%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%p%f(:,:,:)
        f(:,:,n3+1)=b%p%f(:,:,1)
        z(1:n3)=b%p%scales(3)%f(:)
        z(n3+1)=b%p%scales(3)%f(1)+twopi
        deallocate (b%p%f)
        deallocate (b%p%scales(3)%f)
        b%p%dims(3)=n3+1
        b%p%f=>f
        b%p%scales(3)%f=>z
c
        b%nps=b%r%dims(3)
        b%ps=>b%r%scales(3)%f
c
      end if
c
c ****** Set the precision of the B that was read in, based
c ****** on the type of the individual HDF files of the components.
c
      if (verbose.gt.0) then
        write (*,*) 'SET_BTYPE: Setting B precision.'
      end if
c
      if (b%r%hdf32.or.b%t%hdf32.or.b%p%hdf32) then
        b%b_is_32bit=.true.
      else
        b%b_is_32bit=.false.
      end if
c
c ****** Snap the outer (t,p) limits to values that are slightly
c ****** larger than the exact values, if they are close enough.
c ****** This will compensate for the reduced precision inherent in
c ****** magnetic fields read in from 32-bit HDF files, which have
c ****** only ~ 7 digits of accuracy.  In this way, positions
c ****** read from 32-bit HDF files that are near the theta=pi
c ****** and phi=2*pi boundaries will be more likely to end up
c ****** inside the domain.
c
      if (abs(b%ts(b%nts)-pi).lt.eps) then
        if (verbose.gt.0) then
          write (*,*) 'SET_BTYPE: Snapping t.'
        end if
        b%ts(b%nts)=pi+3._r_typ*spacing(pi_r4)
      end if
c
      if (abs(b%ps(b%nps)-twopi).lt.eps) then
        if (verbose.gt.0) then
          write (*,*) 'SET_BTYPE: Snapping p.'
        end if
        b%ps(b%nps)=twopi+3._r_typ*spacing(twopi_r4)
      end if
c
c ****** Set the domain limits.
c
      if (verbose.gt.0) then
        write (*,*) 'SET_BTYPE: Set domain limits.'
      end if
c
      b%lim0(1)=b%rs(1)
      b%lim1(1)=b%rs(b%nrs)
      b%lim0(2)=b%ts(1)
      b%lim1(2)=b%ts(b%nts)
      b%lim0(3)=b%ps(1)
      b%lim1(3)=b%ps(b%nps)
c
c ****** Build the inverse interpolation tables for the
c ****** main mesh.
c
      if (verbose.gt.0) then
        write (*,*) 'SET_BTYPE: Building inverse interpolation tables.'
      end if
c
      b%rs_invtab%n=b%nrs
      allocate (b%rs_invtab%f(b%rs_invtab%n))
      call getinv (b%rs,b%nrs,b%rs_invtab)
c
      b%ts_invtab%n=b%nts
      allocate (b%ts_invtab%f(b%ts_invtab%n))
      call getinv (b%ts,b%nts,b%ts_invtab)
c
      b%ps_invtab%n=b%nps
      allocate (b%ps_invtab%f(b%ps_invtab%n))
      call getinv (b%ps,b%nps,b%ps_invtab)
c
c ****** Compute the mesh cell dimensions on the main mesh.
c
      if (verbose.gt.0) then
        write (*,*) 'SET_BTYPE: Computing mesh cell dims.'
      end if
c
c ****** These are used in setting the field line integration
c ****** step size.
c
      allocate (b%drs(b%nrs))
      allocate (b%dts(b%nts))
      allocate (b%dps(b%nps))
      allocate (b%sts(b%nts))
c
      call get_dx (b%nrs,b%rs,b%drs)
      call get_dx (b%nts,b%ts,b%dts)
      call get_dx (b%nps,b%ps,b%dps)
      b%sts=sin(b%ts)
      b%sts(    1)=max(b%sts(    1),sin(b%dts(    1)))
      b%sts(b%nts)=max(b%sts(b%nts),sin(b%dts(b%nts)))
c
      return
      end
c#######################################################################
      subroutine get_dx (n,x,dx)
c
c-----------------------------------------------------------------------
c
c ****** Get the cell size DX(N) from the 1D mesh in X(N).
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,dx
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      if (n.le.1) then
        dx(1)=0.
      else if (n.eq.2) then
        dx(1)=x(2)-x(1)
        dx(2)=dx(1)
      else
        do i=2,n-1
          dx(i)=half*(x(i+1)-x(i-1))
        enddo
        dx(1)=dx(2)
        dx(n)=dx(n-1)
      end if
c
      return
      end
c#######################################################################
      subroutine make_new_meshes (b)
c
c-----------------------------------------------------------------------
c
c ****** Make new r, t, and p meshes, if requested, or link them
c ****** to the meshes in the B files.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use params
      use files
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      integer :: i,ierr
      real(r_typ) :: d
c
c-----------------------------------------------------------------------
c
c ****** Make the r mesh.
c
      if (new_r_mesh) then
c
c ****** Check if the mesh is to be read from a 1D HDF file.
c
        if (mesh_file_r.ne.' ') then
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Reading the r mesh from file: ',
     &                  trim(mesh_file_r)
          end if
c
          call rdhdf (mesh_file_r,s,ierr)
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the r mesh'//
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_r)
            call exit (1)
          end if
c
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the r mesh'//
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_r)
            call exit (1)
          end if
c
          nrss=s%dims(1)
          allocate (rss(nrss))
          rss=s%f(:,1,1)
c
          call deallocate_sds (s)
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Mesh read in for the r mesh:'
            write (*,*) 'Number of points = ',nrss
            write (*,*) 'Lower limit = ',rss(1)
            write (*,*) 'Upper limit = ',rss(nrss)
          end if
c
        else
c
c ****** Generate a uniform mesh.
c
          if (nrss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//
     &                  ' for the uniform r mesh.'
            write (*,*) 'Number of points specified = ',nrss
            call exit (1)
          end if
c
          allocate (rss(nrss))
c
          if (r0.eq.0..and.r1.eq.0.) then
            r0=b%lim0(1)
            r1=b%lim1(1)
          end if
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Generating a uniform r mesh:'
            write (*,*) 'Number of points = ',nrss
            write (*,*) 'Lower limit = ',r0
            write (*,*) 'Upper limit = ',r1
          end if
c
          if (nrss.ne.1) then
            d=(r1-r0)/(nrss-1)
          else
            d=0.
          end if
c
          do i=1,nrss
            rss(i)=r0+(i-1)*d
          enddo
          rss(1)=r0
          if (nrss.gt.1) then
            rss(nrss)=r1
          end if
c
        end if
c
      else
c
c ****** Use the same mesh as the primary B field mesh.
c
        nrss=b%nrs
        rss=>b%rs
c
      end if
c
c ****** Make the t mesh.
c
      if (new_t_mesh) then
c
c ****** Check if the mesh is to be read from a 1D HDF file.
c
        if (mesh_file_t.ne.' ') then
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Reading the t mesh from file: ',
     &                  trim(mesh_file_t)
          end if
c
          call rdhdf (mesh_file_t,s,ierr)
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
c
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
c
          ntss=s%dims(1)
          allocate (tss(ntss))
          tss=s%f(:,1,1)
c
          call deallocate_sds (s)
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Mesh read in for the t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',tss(1)
            write (*,*) 'Upper limit = ',tss(ntss)
          end if
c
        else
c
c ****** Generate a uniform mesh.
c
          if (ntss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//
     &                  ' for the uniform t mesh.'
            write (*,*) 'Number of points specified = ',ntss
            call exit (1)
          end if
c
          allocate (tss(ntss))
c
          if (t0.eq.0..and.t1.eq.0.) then
            t0=b%lim0(2)
            t1=b%lim1(2)
          end if
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Generating a uniform t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',t0
            write (*,*) 'Upper limit = ',t1
          end if
c
          if (ntss.ne.1) then
            d=(t1-t0)/(ntss-1)
          else
            d=0.
          end if
c
          do i=1,ntss
            tss(i)=t0+(i-1)*d
          enddo
          tss(1)=t0
          if (ntss.gt.1) then
            tss(ntss)=t1
          end if
c
        end if
c
      else
c
c ****** Use the same mesh as the primary B field mesh.
c
        ntss=b%nts
        tss=>b%ts
c
      end if
c
c ****** Make the p mesh.
c
      if (new_p_mesh) then
c
c ****** Check if the mesh is to be read from a 1D HDF file.
c
        if (mesh_file_p.ne.' ') then
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Reading the p mesh from file: ',
     &                  trim(mesh_file_p)
          end if
c
          call rdhdf (mesh_file_p,s,ierr)
c
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
c
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
c
          npss=s%dims(1)
          allocate (pss(npss))
          pss=s%f(:,1,1)
c
          call deallocate_sds (s)
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Mesh read in for the p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',pss(1)
            write (*,*) 'Upper limit = ',pss(npss)
          end if
c
        else
c
c ****** Generate a uniform mesh.
c
          if (npss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//
     &                  ' for the uniform p mesh.'
            write (*,*) 'Number of points specified = ',npss
            call exit (1)
          end if
c
          allocate (pss(npss))
c
          if (p0.eq.0..and.p1.eq.0.) then
            p0=b%lim0(3)
            p1=b%lim1(3)
          end if
c
          if (verbose.gt.0) then
            write (*,*)
            write (*,*) '### Generating a uniform p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',p0
            write (*,*) 'Upper limit = ',p1
          end if
c
          if (npss.ne.1) then
            d=(p1-p0)/(npss-1)
          else
            d=0.
          end if
c
          do i=1,npss
            pss(i)=p0+(i-1)*d
          enddo
          pss(1)=p0
          if (npss.gt.1) then
            pss(npss)=p1
          end if
c
        end if
c
      else
c
c ****** Use the same mesh as the primary B field mesh.
c
        npss=b%nps
        pss=>b%ps
c
      end if
c
      return
      end
c#######################################################################
      subroutine set_ds (b)
c
c-----------------------------------------------------------------------
c
c ****** Set the field line integration step size.
c
c-----------------------------------------------------------------------
c
c ****** If SET_DS_AUTOMATICALLY=.T., the miniumum step size is set
c ****** to the minimum of the cell dimensions from the magnetic
c ****** field files, and the maximum step size is set to the
c ****** maximum of the cell dimensions.  Otherwise, the values read
c ****** in for DS%MIN and DS%MAX are used.
c
c ****** After being set in the above way, DS%MIN and DS%MAX are
c ****** multiplied by the factor DSMULT.  Thus, DSMULT provides a
c ****** quick way to change the integration step size.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use constants
      use vars
      use params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
      real(r_typ) :: dr,dt,dp
      real(r_typ) :: drmin,dtmin,dpmin
      real(r_typ) :: drmax,dtmax,dpmax
c
c-----------------------------------------------------------------------
c
      if (set_ds_automatically) then
c
        drmin=abs(b%lim1(1)-b%lim0(1))
        drmax=0.
        do i=1,b%nrs-1
          dr=abs(b%rs(i+1)-b%rs(i))
          drmin=min(drmin,dr)
          drmax=max(drmax,dr)
        enddo
c
        dtmin=pi
        dtmax=0.
        do j=1,b%nts-1
          dt=abs(b%ts(j+1)-b%ts(j))
          dtmin=min(dtmin,dt)
          dtmax=max(dtmax,dt)
        enddo
c
        dpmin=twopi
        dpmax=0.
        do k=1,b%nps-1
          dp=abs(b%ps(k+1)-b%ps(k))
          dpmin=min(dpmin,dp)
          dpmax=max(dpmax,dp)
        enddo
c
        ds%min=min(drmin,b%lim0(1)*dtmin,b%lim0(1)*dtmin*dpmin)
        ds%max=max(drmax,b%lim1(1)*dtmax,b%lim1(1)*dpmax)
c
      end if
c
      if (dsmult.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in SET_DS:'
        write (*,*) '### DSMULT must be positive.'
        write (*,*) 'DSMULT= ',dsmult
        call exit (1)
      end if
c
      ds%over_rc=ds%over_rc*dsmult
      ds%min=ds%min*dsmult
      ds%max=ds%max*dsmult
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Field line integration parameters:'
        if (ds%variable) then
          write (*,*)
          write (*,*) '### Integration step size control:'//
     &                ' variable step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS%OVER_RC = ',ds%over_rc
          write (*,*) 'DS%MIN = ',ds%min
          write (*,*) 'DS%MAX = ',ds%max
          write (*,*)
          if (ds%limit_by_local_mesh) then
            write (*,*) '### Step size limited by the local'//
     &                  ' B mesh: yes'
            write (*,*) 'DS%LOCAL_MESH_FACTOR = ',
     &                  ds%local_mesh_factor
          else
            write (*,*) '### Step size limited by the local'//
     &                  ' B mesh: no'
          end if
        else
          write (*,*)
          write (*,*) '### Integration step size control:'//
     &                ' uniform step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS = ',ds%min
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine map_forward
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines outward from r=R0.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl,length
c
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl,slogqfl
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines, starting from each (T,P) cell at r=R0,
c ****** until the field line hits r=R1, or goes back to r=R0,
c ****** or exhausts the field line length allowed.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a forward mapping from R0:'
      end if
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      nbad=0
c
      ds%direction_is_along_b=.false.
      ds%direction=1
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose.gt.0) then
c$omp critical
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical
          end if
c
          xfl0(1)=b%lim0(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1, and set
c ****** the expansion factor.
c
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            length(j,k)=s
            if (bs1(1).ne.0.) then
              efl(j,k)=abs((bs0(1)*xfl0(1)**2)/(bs1(1)*xfl1(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
c$omp critical
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_FORWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
c$omp end critical
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            length(j,k)=0.
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
c
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
c$omp critical
              write (*,*)
              write (*,*) '### ERROR in MAP_FORWARD:'
              write (*,*) '### Too many field lines did not reach'//
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
c$omp end critical
            end if
          end if
c
c ****** Write progress diagnostics if requested.
c
          if (verbose.gt.0) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the mapping.
c
      wrote_cr=.false.
c
      if (rffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                'r to file: ',
     &                trim(rffile)
        end if
        call wrhdf_2d (rffile,.true.,ntss,npss,rfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
c
      if (tffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                't to file: ',
     &                trim(tffile)
        end if
        call wrhdf_2d (tffile,.true.,ntss,npss,tfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
c
      if (pffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//
     &                'p to file: ',
     &                trim(pffile)
        end if
        call wrhdf_2d (pffile,.true.,ntss,npss,pfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
c
      if (effile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//
     &                'expansion factor to file: ',
     &                trim(effile)
        end if
        call wrhdf_2d (effile,.true.,ntss,npss,efl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for the expansion factor.'
          call exit (1)
        end if
      end if
c
      if (kffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//
     &                'K factor to file: ',
     &                trim(kffile)
        end if
        call wrhdf_2d (kffile,.true.,ntss,npss,kfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                ' file for the K factor.'
          call exit (1)
        end if
      end if
c
      if (lffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//
     &                'length to file: ',
     &                trim(lffile)
        end if
        call wrhdf_2d (lffile,.true.,ntss,npss,length,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//
     &                'length file.'
          call exit (1)
        end if
      end if
c
c ****** Compute Q (if requested).
c
      if (qffile.eq.' '.and.slogqffile.eq.' ') return
c
c ****** This can only be done if NTSS and NPSS exceed 1.
c
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_FORWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//
     &              ' must be greater than 1.'
        return
      end if
c
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
      if (slogqffile.ne.' ') then
        allocate (slogqfl(ntss-1,npss-1))
      end if
c
c ****** Define the half-meshes (on which Q is computed).
c
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
c
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
c
ccc$omp parallel do
ccc$omp& default(private)
ccc$omp& shared(efl,pfl,tfl,qfl,pss,tss,tssh)
ccc$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.
     &        efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)/efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
c
      if (qffile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//
     &                'Q to file: ',
     &                trim(qffile)
        end if
        call wrhdf_2d (qffile,.true.,ntss-1,npss-1,qfl,tssh,pssh,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the Q factor file.'
          call exit (1)
        end if
      end if
c
      if (slogqffile.ne.' ') then
c
        call slogq (qfl,slogqfl,b%lim0(1))
c
        if (verbose.gt.0) then
          write (*,*) 'Writing the forward mapping '//
     &                'SLOG(Q) to file: ',
     &                trim(slogqffile)
        end if
        call wrhdf_2d(slogqffile,.true.,ntss-1,npss-1,slogqfl,tssh,pssh,
     &                hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the SLOG(Q) factor file.'
          call exit (1)
        end if
        deallocate (slogqfl)
      end if
c
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
c
      return
      end
c#######################################################################
      subroutine slogq (qfl,slogqfl,rlevel)
c
c-----------------------------------------------------------------------
c
c ****** Compute Slava's "signed log" of Q.
c
c-----------------------------------------------------------------------
c
      use number_types
      use field
      use mesh
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1.0_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(ntss-1,npss-1) :: qfl,slogqfl
      real(r_typ) :: rlevel
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(ntss-1) :: tssh
      real(r_typ), dimension(npss-1) :: pssh
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
      real(r_typ) :: q,lq,br
      type(csvec) :: x,bv
c
c-----------------------------------------------------------------------
c
c ****** Define the half-meshes (on which Q was computed).
c
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
c
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
c
c ***** Calculate slogq
c
!$omp parallel do collapse(2)
!$omp& default(shared) private(i,j,q,lq,x,br,bv)
      do i=1,ntss-1
        do j=1,npss-1
          q=half*qfl(i,j)
          if (q.lt.one) then
            q=one
          else
            q=q+sqrt(q**2-one)
          end if
c
          lq=log10(q)
c
c ****** Get Br at the current location.
c
          x%s(1)=rlevel
          x%s(2)=tssh(i)
          x%s(3)=pssh(j)
          call getb (b,x,bv)
          br=bv%s(1)
          if (br.lt.0) then
            lq=-lq
          end if
          slogqfl(i,j)=lq
        enddo
      enddo
!$omp end parallel do
c
      end subroutine
c#######################################################################
      subroutine map_backward
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines inward from r=R1.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl,length
c
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl,slogqfl
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines, starting from each (T,P) cell at r=R1,
c ****** until the field line hits r=R0, or goes back to r=R1,
c ****** or exhausts the number of segments allowed.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a backward mapping from R1:'
      end if
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      nbad=0
c
      ds%direction_is_along_b=.false.
      ds%direction=-1
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& default(shared)
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose.gt.0) then
c$omp critical (omp_nc)
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical (omp_nc)
          end if
c
          xfl0(1)=b%lim1(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
c
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            length(j,k)=s
            if (bs0(1).ne.0.) then
              efl(j,k)=abs((bs1(1)*xfl1(1)**2)/(bs0(1)*xfl0(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
c$omp critical (nbad_count)
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_BACKWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
c$omp end critical (nbad_count)
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            length(j,k)=0.
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
c
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
c$omp critical (nbad2)
              write (*,*)
              write (*,*) '### ERROR in MAP_BACKWARD:'
              write (*,*) '### Too many field lines did not reach'//
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
c$omp end critical (nbad2)
            end if
          end if
c
c ****** Write progress diagnostics if requested.
c
          if (verbose.gt.0) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the mapping.
c
      wrote_cr=.false.
c
      if (rbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                'r to file: ',
     &                trim(rbfile)
        end if
        call wrhdf_2d (rbfile,.true.,ntss,npss,rfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
c
      if (tbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                't to file: ',
     &                trim(tbfile)
        end if
        call wrhdf_2d (tbfile,.true.,ntss,npss,tfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
c
      if (pbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//
     &                'p to file: ',
     &                trim(pbfile)
        end if
        call wrhdf_2d (pbfile,.true.,ntss,npss,pfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
c
      if (ebfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//
     &                'expansion factor to file: ',
     &                trim(ebfile)
        end if
        call wrhdf_2d (ebfile,.true.,ntss,npss,efl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for the expansion factor.'
          call exit (1)
        end if
      end if
c
      if (kbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//
     &                'K factor to file: ',
     &                trim(kbfile)
        end if
        call wrhdf_2d (kbfile,.true.,ntss,npss,kfl,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                ' file for the K factor.'
          call exit (1)
        end if
      end if
c
      if (lbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//
     &                'length to file: ',
     &                trim(lbfile)
        end if
        call wrhdf_2d (lbfile,.true.,ntss,npss,length,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//
     &                'length file.'
          call exit (1)
        end if
      end if
c
c ****** Compute Q (if requested).
c
      if (qbfile.eq.' '.and.slogqbfile.eq.' ') return
c
c ****** This can only be done if NTSS and NPSS exceed 1.
c
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_BACKWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//
     &              ' must be greater than 1.'
        return
      end if
c
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
      if (slogqbfile.ne.' ') then
        allocate (slogqfl(ntss-1,npss-1))
      end if
c
c ****** Define the half-meshes (on which Q is computed).
c
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
c
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
c
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.
     &        efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)*efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
c
      if (qbfile.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//
     &                'Q to file: ',
     &                trim(qbfile)
        end if
        call wrhdf_2d (qbfile,.true.,ntss-1,npss-1,qfl,tssh,pssh,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the Q factor file.'
          call exit (1)
        end if
      end if
c
      if (slogqbfile.ne.' ') then
c
        call slogq (qfl,slogqfl,b%lim1(1))
c
        if (verbose.gt.0) then
          write (*,*) 'Writing the backward mapping '//
     &                'SLOG(Q) to file: ',
     &                trim(slogqbfile)
        end if
        call wrhdf_2d(slogqbfile,.true.,ntss-1,npss-1,slogqfl,tssh,pssh,
     &                hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the SLOG(Q) factor file.'
          call exit (1)
        end if
        deallocate (slogqfl)
      end if
c
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
c
      return
      end
c#######################################################################
      function modulo_twopi (x)
c
c-----------------------------------------------------------------------
c
c ****** Return the smallest value of X, modulo 2*pi.
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x
      real(r_typ) :: modulo_twopi
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xm,xp,x_min
c
c-----------------------------------------------------------------------
c
      xm=abs(x-twopi)
      xp=abs(x+twopi)
      x_min=min(xm,abs(x),xp)
      if (xm.eq.x_min) then
        modulo_twopi=x-twopi
      else if (xp.eq.x_min) then
        modulo_twopi=x+twopi
      else
        modulo_twopi=x
      end if
c
      return
      end
c#######################################################################
      subroutine map_3d
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines from every point on a 3D (r,t,p) mesh.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(nrss,ntss,npss) :: rfl,tfl,pfl
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines, starting from each (R,T,P) cell,
c ****** until the field line hits the boundaries or exhausts
c ****** the field line length allowed.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a mapping in 3D:'
      end if
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      ds%direction_is_along_b=.false.
      ds%direction=-1
c
      n_total=nrss*ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(i,j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
          do i=1,nrss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose.gt.0) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
c
            call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
c
            if (ttb) then
              rfl(i,j,k)=xfl1(1)
            else
              rfl(i,j,k)=-xfl1(1)
            end if
            tfl(i,j,k)=xfl1(2)
            pfl(i,j,k)=xfl1(3)
c
c ****** Write progress diagnostics if requested.
c
            if (verbose.gt.0) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the mapping.
c
      wrote_cr=.false.
c
      if (volume3d_output_file%r.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//
     &                'r to file: ',
     &                trim(volume3d_output_file%r)
        end if
        call wrhdf_3d (volume3d_output_file%r,.true.,
     &                 nrss,ntss,npss,rfl,rss,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
c
      if (volume3d_output_file%t.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//
     &                't to file: ',
     &                trim(volume3d_output_file%t)
        end if
        call wrhdf_3d (volume3d_output_file%t,.true.,
     &                 nrss,ntss,npss,tfl,rss,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
c
      if (volume3d_output_file%p.ne.' ') then
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//
     &                'p to file: ',
     &                trim(volume3d_output_file%p)
        end if
        call wrhdf_3d (volume3d_output_file%p,.true.,
     &                 nrss,ntss,npss,pfl,rss,tss,pss,
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine read_slice_coordinates
c
c-----------------------------------------------------------------------
c
c ****** Read the coordinates that define the slice in the 3D volume.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      logical, external :: same_structure_sds
c
c-----------------------------------------------------------------------
c
c ****** Set the coordinate names.
c
      if (slice_coords_are_xyz) then
        slice_coord_name=(/'x','y','z'/)
      else
        slice_coord_name=(/'r','t','p'/)
      end if
c
c ****** Read the coordinates of the slice.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Reading the coordinates of the slice ...'
      end if
c
c ****** Read the x/r coordinate file.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Reading the '//slice_coord_name(1)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%r)
      end if
c
      call rdhdf (slice_input_file%r,slice_c1,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(1)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%r)
        call exit (1)
      end if
c
c ****** Read the y/t coordinate file.
c
      if (verbose.gt.0) then
        write (*,*) 'Reading the '//slice_coord_name(2)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%t)
      end if
c
      call rdhdf (slice_input_file%t,slice_c2,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(2)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%t)
        call exit (1)
      end if
c
c ****** Check that the y/t coordinate has the same structure as the
c ****** x/r coordinate.
c
      if (.not.same_structure_sds(slice_c1,slice_c2)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//
     &              slice_coord_name(1)//' and '//
     &              slice_coord_name(2)//' do not have'//
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(2)//
     &              ' file name: ',trim(slice_input_file%t)
        call exit (1)
      end if
c
c ****** Read the z/p coordinate file.
c
      if (verbose.gt.0) then
        write (*,*) 'Reading the '//slice_coord_name(3)//
     &              ' coordinate from file: ',
     &              trim(slice_input_file%p)
      end if
c
      call rdhdf (slice_input_file%p,slice_c3,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(3)//
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%p)
        call exit (1)
      end if
c
c ****** Check that the z/p coordinate has the same structure as the
c ****** x/r coordinate.
c
      if (.not.same_structure_sds(slice_c1,slice_c3)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//
     &              slice_coord_name(1)//' and '//
     &              slice_coord_name(3)//' do not have'//
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(3)//
     &              ' file name: ',trim(slice_input_file%p)
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine deallocate_slice_coordinates
c
c-----------------------------------------------------------------------
c
      use sds_def
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Read the coordinates of the slice.
c
      call deallocate_sds (slice_c1)
      call deallocate_sds (slice_c2)
      call deallocate_sds (slice_c3)
c
      return
      end
c#######################################################################
      function same_structure_sds (s1,s2)
c
c-----------------------------------------------------------------------
c
c ****** Check if the two data sets S1 and S2 have the same
c ****** structure.  If they do, return .TRUE; otherwise, return
c ****** .FALSE.
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s1,s2
      logical :: same_structure_sds
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      same_structure_sds=.false.
c
      if (s1%ndim.ne.s2%ndim) return
c
      if (s1%scale.neqv.s2%scale) return
c
      do i=1,s1%ndim
        if (s1%dims(i).ne.s2%dims(i)) return
      enddo
c
      same_structure_sds=.true.
c
      return
      end
c#######################################################################
      subroutine map_slice
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines from points on a slice in the 3D volume.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
      use diags
      use openmp_vars
      use tracefl_interface
      use debug
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage for the forward mapping.
c
      real(r_typ), dimension(:,:,:), allocatable, target :: rfl_f
      real(r_typ), dimension(:,:,:), allocatable, target :: tfl_f
      real(r_typ), dimension(:,:,:), allocatable, target :: pfl_f
c
c ****** Storage for the backward mapping.
c
      real(r_typ), dimension(:,:,:), allocatable, target :: rfl_b
      real(r_typ), dimension(:,:,:), allocatable, target :: tfl_b
      real(r_typ), dimension(:,:,:), allocatable, target :: pfl_b
c
c-----------------------------------------------------------------------
c
      type (sds) :: out
      integer :: ierr,i,j,k,n1,n2,n3
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      real(r_typ), dimension(3) :: c
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
c ****** Field line trace storage buffers.
c
      type(traj), dimension(:,:,:), allocatable :: xtf
      type(traj), dimension(:,:,:), allocatable :: xtb
      real(r_typ), dimension(:,:), allocatable :: fl
      real(r_typ), dimension(3) :: xyz
      integer :: n,l
      real(r_typ) :: dummy
      character(4) :: ch4
      character(256) :: fname
c
c-----------------------------------------------------------------------
c
c ****** Map the field lines for all points on a slice in 3D.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a mapping starting on a slice:'
      end if
c
c ****** Check that at least one tracing direction was requested.
c
      if (.not.(trace_from_slice_forward.or.
     &          trace_from_slice_backward)) then
        write (*,*)
        write (*,*) '### ERROR in MAP_SLICE:'
        write (*,*) '### At least one tracing direction (forward'//
     &              ' and/or backward) must be requested.'
        call exit (1)
      end if
c
c ****** Set the tracing direction to be either along the direction
c ****** of the magnetic field or along the directon of increasing
c ****** radius.
c
      ds%direction_is_along_b=trace_slice_direction_is_along_b
c
      ds_f=ds
      ds_f%direction=1
c
      ds_b=ds
      ds_b%direction=-1
c
      if (verbose.gt.0) then
        if (trace_slice_direction_is_along_b) then
          write (*,*)
          write (*,*) '### The forward tracing direction is'//
     &                ' along the direction of B.'
          write (*,*) '### The backward tracing direction is'//
     &                ' opposite to the direction of B.'
        else
          write (*,*)
          write (*,*) '### The forward tracing direction is'//
     &                ' along the direction of increasing radius.'
          write (*,*) '### The backward tracing direction is'//
     &                ' along the direction of decreasing radius.'
        end if
      end if
c
      n1=slice_c1%dims(1)
      n2=slice_c1%dims(2)
      n3=slice_c1%dims(3)
c
c ****** Allocate storage for the mapping.
c
      if (trace_from_slice_forward) then
        allocate (rfl_f(n1,n2,n3))
        allocate (tfl_f(n1,n2,n3))
        allocate (pfl_f(n1,n2,n3))
      end if
c
      if (trace_from_slice_backward) then
        allocate (rfl_b(n1,n2,n3))
        allocate (tfl_b(n1,n2,n3))
        allocate (pfl_b(n1,n2,n3))
      end if
c
c ****** Allocate the field line storage buffers if the field
c ****** line traces are being written to HDF files.
c
      if (write_traces_to_hdf) then
        if (trace_from_slice_forward) then
          allocate (xtf(n1,n2,n3))
          do k=1,n3
            do j=1,n2
              do i=1,n1
                call allocate_trajectory_buffer (xtf(i,j,k))
              enddo
            enddo
          enddo
        end if
        if (trace_from_slice_backward) then
          allocate (xtb(n1,n2,n3))
          do k=1,n3
            do j=1,n2
              do i=1,n1
                call allocate_trajectory_buffer (xtb(i,j,k))
              enddo
            enddo
          enddo
        end if
      end if
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
c ****** Trace a field line from each point on the slice until the
c ****** field line hits the boundaries or exhausts the length
c ****** allowed.
c
      n_total=n1*n2*n3
      n_completed=0
c
c$omp parallel do
c$omp& private(i,j,k,c,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,n3
        do j=1,n2
          do i=1,n1
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose.gt.0) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            if (slice_coords_are_xyz) then
              c=(/slice_c1%f(i,j,k),
     &            slice_c2%f(i,j,k),
     &            slice_c3%f(i,j,k)/)
              call c2s (c,xfl0)
            else
              xfl0=(/slice_c1%f(i,j,k),
     &               slice_c2%f(i,j,k),
     &               slice_c3%f(i,j,k)/)
            end if
c
c ****** Launch a field line in the positive direction, if requested.
c
            if (trace_from_slice_forward) then
c
              if (write_traces_to_hdf) then
                call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb,
     &                        xtf(i,j,k))
              else
                call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
              end if
c
              if (slice_coords_are_xyz) then
                if (ttb) then
                  call s2c (xfl1,c)
                  rfl_f(i,j,k)=c(1)
                  tfl_f(i,j,k)=c(2)
                  pfl_f(i,j,k)=c(3)
                else
                  rfl_f(i,j,k)=0.
                  tfl_f(i,j,k)=0.
                  pfl_f(i,j,k)=0.
                end if
              else
                if (ttb) then
                  rfl_f(i,j,k)=xfl1(1)
                else
                  rfl_f(i,j,k)=-xfl1(1)
                end if
                tfl_f(i,j,k)=xfl1(2)
                pfl_f(i,j,k)=xfl1(3)
              end if
c
            end if
c
c ****** Launch a field line in the negative direction, if requested.
c
            if (trace_from_slice_backward) then
c
              if (write_traces_to_hdf) then
                call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb,
     &                        xtb(i,j,k))
              else
                call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
              end if
c
              if (slice_coords_are_xyz) then
                if (ttb) then
                  call s2c (xfl1,c)
                  rfl_b(i,j,k)=c(1)
                  tfl_b(i,j,k)=c(2)
                  pfl_b(i,j,k)=c(3)
                else
                  rfl_b(i,j,k)=0.
                  tfl_b(i,j,k)=0.
                  pfl_b(i,j,k)=0.
                end if
              else
                if (ttb) then
                  rfl_b(i,j,k)=xfl1(1)
                else
                  rfl_b(i,j,k)=-xfl1(1)
                end if
                tfl_b(i,j,k)=xfl1(2)
                pfl_b(i,j,k)=xfl1(3)
              end if
c
            end if
c
c ****** Write progress diagnostics if requested.
c
            if (verbose.gt.0) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the forward mapping.
c
      wrote_cr=.false.
c
      if (trace_from_slice_forward.and.
     &    slice_output_file_forward%r.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>rfl_f
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//
     &                ' coordinate '//slice_coord_name(1)//
     &                ' to file: ',
     &                trim(slice_output_file_forward%r)
        end if
        call wrhdf (slice_output_file_forward%r,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(1)//'.'
          call exit (1)
        end if
      end if
c
      if (trace_from_slice_forward.and.
     &    slice_output_file_forward%t.ne.' ') then
        out%ndim=slice_c2%ndim
        out%dims=slice_c2%dims
        out%scale=slice_c2%scale
        out%hdf32=slice_c2%hdf32
        out%scales(1)%f=>slice_c2%scales(1)%f
        out%scales(2)%f=>slice_c2%scales(2)%f
        out%scales(3)%f=>slice_c2%scales(3)%f
        out%f=>tfl_f
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//
     &                ' coordinate '//slice_coord_name(2)//
     &                ' to file: ',
     &                trim(slice_output_file_forward%t)
        end if
        call wrhdf (slice_output_file_forward%t,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(2)//'.'
          call exit (1)
        end if
      end if
c
      if (trace_from_slice_forward.and.
     &    slice_output_file_forward%p.ne.' ') then
        out%ndim=slice_c3%ndim
        out%dims=slice_c3%dims
        out%scale=slice_c3%scale
        out%hdf32=slice_c3%hdf32
        out%scales(1)%f=>slice_c3%scales(1)%f
        out%scales(2)%f=>slice_c3%scales(2)%f
        out%scales(3)%f=>slice_c3%scales(3)%f
        out%f=>pfl_f
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//
     &                ' coordinate '//slice_coord_name(3)//
     &                ' to file: ',
     &                trim(slice_output_file_forward%p)
        end if
        call wrhdf (slice_output_file_forward%p,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(3)//'.'
          call exit (1)
        end if
      end if
c
c ****** Write the backward mapping.
c
      wrote_cr=.false.
c
      if (trace_from_slice_backward.and.
     &    slice_output_file_backward%r.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>rfl_b
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//
     &                ' coordinate '//slice_coord_name(1)//
     &                ' to file: ',
     &                trim(slice_output_file_backward%r)
        end if
        call wrhdf (slice_output_file_backward%r,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(1)//'.'
          call exit (1)
        end if
      end if
c
      if (trace_from_slice_backward.and.
     &    slice_output_file_backward%t.ne.' ') then
        out%ndim=slice_c2%ndim
        out%dims=slice_c2%dims
        out%scale=slice_c2%scale
        out%hdf32=slice_c2%hdf32
        out%scales(1)%f=>slice_c2%scales(1)%f
        out%scales(2)%f=>slice_c2%scales(2)%f
        out%scales(3)%f=>slice_c2%scales(3)%f
        out%f=>tfl_b
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//
     &                ' coordinate '//slice_coord_name(2)//
     &                ' to file: ',
     &                trim(slice_output_file_backward%t)
        end if
        call wrhdf (slice_output_file_backward%t,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(2)//'.'
          call exit (1)
        end if
      end if
c
      if (trace_from_slice_backward.and.
     &    slice_output_file_backward%p.ne.' ') then
        out%ndim=slice_c3%ndim
        out%dims=slice_c3%dims
        out%scale=slice_c3%scale
        out%hdf32=slice_c3%hdf32
        out%scales(1)%f=>slice_c3%scales(1)%f
        out%scales(2)%f=>slice_c3%scales(2)%f
        out%scales(3)%f=>slice_c3%scales(3)%f
        out%f=>pfl_b
        if (verbose.gt.0) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//
     &                ' coordinate '//slice_coord_name(3)//
     &                ' to file: ',
     &                trim(slice_output_file_backward%p)
        end if
        call wrhdf (slice_output_file_backward%p,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//
     &                ' mapping file for coordinate '//
     &                slice_coord_name(3)//'.'
          call exit (1)
        end if
      end if
c
c ****** Write the forward field line traces to individual HDF
c ****** files if requested.
c
      if (write_traces_to_hdf.and.trace_from_slice_forward) then
c
        if (verbose.gt.0) then
          write (*,*)
        end if
c
c ****** Loop over all points.
c
        n=0
        do k=1,n3
          do j=1,n2
            do i=1,n1
c
c ****** Allocate a temporary array to store the field line trace.
c
              allocate (fl(3,xtf(i,j,k)%npts))
c
c ****** Load the array with the field line coordinates.
c
              do l=1,xtf(i,j,k)%npts
                fl(1,l)=xtf(i,j,k)%x(1)%f(l)
                fl(2,l)=xtf(i,j,k)%x(2)%f(l)
                fl(3,l)=xtf(i,j,k)%x(3)%f(l)
                if (write_traces_as_xyz) then
                  call s2c (fl(1,l),xyz)
                  fl(1:3,l)=xyz
                end if
              enddo
c
c ****** Write the HDF files.
c
              n=n+1
              write (ch4,'(i4.4)') n
              fname=trim(write_traces_root)//'_f_'//ch4//'.'//trim(fmt)
c
              if (verbose.gt.0) then
                write (*,*) 'Writing a forward field line trace'//
     &                      ' to file: ',trim(fname)
              end if
c
              call wrhdf_2d (fname,.false.,
     &                       3,xtf(i,j,k)%npts,fl,
     &                       fl(:,1),fl(1,:),
     &                       slice_c1%hdf32,ierr)
c
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ERROR in MAP_SLICE:'
                write (*,*) '### Could not write a forward'//
     &                      ' field line trace to HDF file:'
                write (*,*) trim(fname)
                call exit (1)
              end if
c
              deallocate (fl)
c
              call deallocate_trajectory_buffer (xtf(i,j,k))
c
            enddo
          enddo
        enddo
c
        deallocate (xtf)
c
      end if
c
c ****** Write the backward field line traces to individual HDF
c ****** files if requested.
c
      if (write_traces_to_hdf.and.trace_from_slice_backward) then
c
        if (verbose.gt.0) then
          write (*,*)
        end if
c
c ****** Loop over all points.
c
        n=0
        do k=1,n3
          do j=1,n2
            do i=1,n1
c
c ****** Allocate a temporary array to store the field line trace.
c
              allocate (fl(3,xtb(i,j,k)%npts))
c
c ****** Load the array with the field line coordinates.
c
              do l=1,xtb(i,j,k)%npts
                fl(1,l)=xtb(i,j,k)%x(1)%f(l)
                fl(2,l)=xtb(i,j,k)%x(2)%f(l)
                fl(3,l)=xtb(i,j,k)%x(3)%f(l)
                if (write_traces_as_xyz) then
                  call s2c (fl(1,l),xyz)
                  fl(1:3,l)=xyz
                end if
              enddo
c
c ****** Write the HDF files.
c
              n=n+1
              write (ch4,'(i4.4)') n
              fname=trim(write_traces_root)//'_b_'//ch4//'.'//trim(fmt)
c
              if (verbose.gt.0) then
                write (*,*) 'Writing a backward field line trace'//
     &                      ' to file: ',trim(fname)
              end if
c
              call wrhdf_2d (fname,.false.,
     &                       3,xtb(i,j,k)%npts,fl,
     &                       fl(:,1),fl(1,:),
     &                       slice_c1%hdf32,ierr)
c
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ERROR in MAP_SLICE:'
                write (*,*) '### Could not write a backward'//
     &                      ' field line trace to HDF file:'
                write (*,*) trim(fname)
                call exit (1)
              end if
c
              deallocate (fl)
c
              call deallocate_trajectory_buffer (xtb(i,j,k))
c
            enddo
          enddo
        enddo
c
        deallocate (xtb)
c
      end if
c
c ****** Deallocate memory.
c
      if (trace_from_slice_forward) then
        deallocate (rfl_f)
        deallocate (tfl_f)
        deallocate (pfl_f)
      end if
c
      if (trace_from_slice_backward) then
        deallocate (rfl_b)
        deallocate (tfl_b)
        deallocate (pfl_b)
      end if
c
      return
      end
c#######################################################################
      subroutine get_ch_map (rv)
c
c-----------------------------------------------------------------------
c
c ****** Compute a coronal hole map at radius r=RV.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the coronal hole map.
c
      real(r_typ), dimension(npss,ntss) :: ch
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a coronal hole map at r = ',rv
      end if
c
c ****** Check that the radius specified is valid.
c
      if (rv.lt.b%lim0(1).or.rv.gt.b%lim1(1)) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Invalid radius specified.'
        write (*,*) '### The radius is outside the domain limits:'
        write (*,*) 'Lower radial domain limit = ',b%lim0(1)
        write (*,*) 'Upper radial domain limit = ',b%lim1(1)
        write (*,*) 'Specified radius          = ',rv
        call exit (1)
      end if
c
c ****** Check that the coronal hole map output file name is not
c ****** blank, since this does not make sense.
c
      if (ch_map_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### A coronal hole map was requested, yet'//
     &              ' the output file name is blank.'
        call exit (1)
      end if
c
c ****** Set the tracing direction to be along the direction
c ****** of the magnetic field.
c
      ds%direction_is_along_b=.true.
c
      ds_f=ds
      ds_f%direction=1
c
      ds_b=ds
      ds_b%direction=-1
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      n_total=ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(f_trace_reached_boundary,f_br_positive)
c$omp& private(f_trace_on_r0,f_trace_on_r1)
c$omp& private(b_trace_reached_boundary,b_br_positive)
c$omp& private(b_trace_on_r0,b_trace_on_r1)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
          if (verbose.gt.0) then
c$omp critical
            n_completed=n_completed+1
            nc=n_completed
c$omp end critical
          end if
c
          xfl0(1)=rv
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
c
c ****** Trace a field line in the forward direction along B.
c
          call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
          if (ttb) then
            f_trace_reached_boundary=.true.
            f_trace_on_r0=xfl1(1).eq.b%lim0(1)
            f_trace_on_r1=xfl1(1).eq.b%lim1(1)
            f_br_positive=bs1(1).ge.0.
          else
            f_trace_reached_boundary=.false.
          end if
c
c ****** Trace a field line in the backward direction along B.
c
          call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
          if (ttb) then
            b_trace_reached_boundary=.true.
            b_trace_on_r0=xfl1(1).eq.b%lim0(1)
            b_trace_on_r1=xfl1(1).eq.b%lim1(1)
            b_br_positive=bs1(1).ge.0.
          else
            b_trace_reached_boundary=.false.
          end if
c
c ****** Set the coronal hole map value.
c
c ****** Note that the following values are set in the output
c ****** coronal hole map:
c ******
c ******    -1: open field line with negative polarity
c ******     1: open field line with positive polarity
c ******     0: closed field line with both footpoints
c ******        on r=R0
c ******     2: closed field line with both footpoints
c ******        on r=R1
c ******    -2: field line that does not reach either
c ******        the r=R0 boundary or the r=R1 boundary
c
          if (f_trace_reached_boundary.and.
     &        b_trace_reached_boundary) then
            if (f_trace_on_r0.and.b_trace_on_r1) then
              if (f_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r1.and.b_trace_on_r0) then
              if (b_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r0.and.b_trace_on_r0) then
              ch(k,j)=0.
            else if (f_trace_on_r1.and.b_trace_on_r1) then
              ch(k,j)=two
            else
              ch(k,j)=-two
            end if
          else
            ch(k,j)=-two
          end if
c
c ****** Write progress diagnostics if requested.
c
          if (verbose.gt.0) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
c
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the coronal hole map.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Writing the coronal hole map to file: ',
     &              trim(ch_map_output_file)
      end if
c
      call wrhdf_2d (ch_map_output_file,.true.,
     &               npss,ntss,ch,pss,tss,
     &               hdf32,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine get_ch_map_3d
c
c-----------------------------------------------------------------------
c
c ****** Compute a 3D coronal hole map.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Storage for the coronal hole map.
c
      real(r_typ), dimension(nrss,ntss,npss) :: ch
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
c
c-----------------------------------------------------------------------
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a 3D coronal hole map:'
      end if
c
c ****** Check that the coronal hole map output file name is not
c ****** blank, since this does not make sense.
c
      if (ch_map_3d_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP_3D:'
        write (*,*) '### A coronal hole map was requested, yet'//
     &              ' the output file name is blank.'
        call exit (1)
      end if
c
c ****** Set the tracing direction to be along the direction
c ****** of the magnetic field.
c
      ds%direction_is_along_b=.true.
c
      ds_f=ds
      ds_f%direction=1
c
      ds_b=ds
      ds_b%direction=-1
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      n_total=nrss*ntss*npss
      n_completed=0
c
c$omp parallel do
c$omp& private(i,j,k,xfl0,xfl1,bs0,bs1,s,ttb)
c$omp& private(f_trace_reached_boundary,f_br_positive)
c$omp& private(f_trace_on_r0,f_trace_on_r1)
c$omp& private(b_trace_reached_boundary,b_br_positive)
c$omp& private(b_trace_on_r0,b_trace_on_r1)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do i=1,nrss
        do k=1,npss
          do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose.gt.0) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
c
c ****** Trace a field line in the forward direction along B.
c
            call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
            if (ttb) then
              f_trace_reached_boundary=.true.
              f_trace_on_r0=xfl1(1).eq.b%lim0(1)
              f_trace_on_r1=xfl1(1).eq.b%lim1(1)
              f_br_positive=bs1(1).ge.0.
            else
              f_trace_reached_boundary=.false.
            end if
c
c ****** Trace a field line in the backward direction along B.
c
            call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
c
c ****** Check that the field line reached R0 or R1.
c
            if (ttb) then
              b_trace_reached_boundary=.true.
              b_trace_on_r0=xfl1(1).eq.b%lim0(1)
              b_trace_on_r1=xfl1(1).eq.b%lim1(1)
              b_br_positive=bs1(1).ge.0.
            else
              b_trace_reached_boundary=.false.
            end if
c
c ****** Set the coronal hole map value.
c
c ****** Note that the following values are set in the output
c ****** coronal hole map:
c ******
c ******    -1: open field line with negative polarity
c ******     1: open field line with positive polarity
c ******     0: closed field line with both footpoints
c ******        on r=R0
c ******     2: closed field line with both footpoints
c ******        on r=R1
c ******    -2: field line that does not reach either
c ******        the r=R0 boundary or the r=R1 boundary
c
            if (f_trace_reached_boundary.and.
     &        b_trace_reached_boundary) then
              if (f_trace_on_r0.and.b_trace_on_r1) then
                if (f_br_positive) then
                  ch(i,j,k)=one
                else
                  ch(i,j,k)=-one
                end if
              else if (f_trace_on_r1.and.b_trace_on_r0) then
                if (b_br_positive) then
                  ch(i,j,k)=one
                else
                  ch(i,j,k)=-one
                end if
              else if (f_trace_on_r0.and.b_trace_on_r0) then
                ch(i,j,k)=0.
              else if (f_trace_on_r1.and.b_trace_on_r1) then
                ch(i,j,k)=two
              else
                ch(i,j,k)=-two
              end if
            else
              ch(i,j,k)=-two
            end if
c
c ****** Write progress diagnostics if requested.
c
            if (verbose.gt.0) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do
c
c ****** Write the coronal hole map.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Writing the 3D coronal hole map to file: ',
     &              trim(ch_map_3d_output_file)
      end if
c
      call wrhdf_3d (ch_map_3d_output_file,.true.,
     &               nrss,ntss,npss,ch,rss,tss,pss,
     &               hdf32,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP_3D:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine get_q_on_slice
c
c-----------------------------------------------------------------------
c
c ****** Trace field lines from points on a slice in the 3D volume,
c ****** getting Q at each point.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
      use diags
      use openmp_vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage for the mapping.
c
      real(r_typ), dimension(:,:,:), allocatable, target :: qfl
      real(r_typ), dimension(:,:,:), allocatable, target :: length
c
c-----------------------------------------------------------------------
c
      type (sds) :: out
      integer :: ierr,i,j,k,n1,n2,n3
      real(r_typ), dimension(3) :: s0
      real(r_typ), dimension(3) :: c
      real(r_typ) :: q
      logical :: gotq
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
      logical :: save_field_line_length
      real(r_typ) :: lfl
      integer :: iseq
c
c-----------------------------------------------------------------------
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing Q on a slice:'
      end if
c
      n1=slice_c1%dims(1)
      n2=slice_c1%dims(2)
      n3=slice_c1%dims(3)
c
c ****** Set a flag to indicate if the field line length
c ****** was requested.
c
      if (slice_length_output_file.ne.' ') then
        save_field_line_length=.true.
      else
        save_field_line_length=.false.
      end if
c
c ****** Allocate the storage for the output Q array.
c
      allocate (qfl(n1,n2,n3))
c
c ****** Allocate the storage for the output field line
c ****** length, if requested.
c
      if (save_field_line_length) then
        allocate (length(n1,n2,n3))
      end if
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Getting Q ...'
        write (*,*)
      end if
c
c ****** Calculate Q at each point on the slice.
c
      ds%direction_is_along_b=.true.
      ds%direction=1
c
      n_total=n1*n2*n3
      n_completed=0
c
c$omp parallel do
c$omp& private(i,j,k,iseq,c,s0,q,gotq,lfl)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do k=1,n3
        do j=1,n2
          do i=1,n1
c
c ****** Compute a sequence number that is used to construct
c ****** file names when diagnostic field line traces are
c ****** being output to files.
c
            iseq=i+(j-1)*n1+(k-1)*n1*n2
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose.gt.0) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            if (slice_coords_are_xyz) then
              c=(/slice_c1%f(i,j,k),
     &            slice_c2%f(i,j,k),
     &            slice_c3%f(i,j,k)/)
              call c2s (c,s0)
            else
              s0=(/slice_c1%f(i,j,k),
     &             slice_c2%f(i,j,k),
     &             slice_c3%f(i,j,k)/)
            end if
c
            call getq (iseq,ds,s0,q_increment_h,q,gotq,lfl)
c
            if (gotq) then
              qfl(i,j,k)=q
            else
              qfl(i,j,k)=0.
            end if
c
            if (save_field_line_length) then
              length(i,j,k)=lfl
            end if
c
c ****** Write progress diagnostics if requested.
c
            if (verbose.gt.0) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do
c
      if (verbose.gt.0) then
        write (*,*)
      end if
c
c ****** Write the Q slice.
c
      if (slice_q_output_file.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>qfl
        if (verbose.gt.0) then
          write (*,*) 'Writing Q in the slice to file: ',
     &                trim(slice_q_output_file)
        end if
        call wrhdf (slice_q_output_file,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in GET_Q_ON_SLICE:'
          write (*,*) '### Could not write Q in the slice'//
     &                ' to file: ',trim(slice_q_output_file)
          call exit (1)
        end if
      end if
c
      deallocate (qfl)
c
c ****** Write the field line length, if requested.
c
      if (save_field_line_length) then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>length
        if (verbose.gt.0) then
          write (*,*) 'Writing field line length in'//
     &                ' the slice to file: ',
     &                trim(slice_length_output_file)
        end if
        call wrhdf (slice_length_output_file,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in GET_Q_ON_SLICE:'
          write (*,*) '### Could not write the field line'//
     &                ' length in the slice to file: ',
     &                trim(slice_length_output_file)
          call exit (1)
        end if
c
        deallocate (length)
c
      end if
c
      return
      end
c#######################################################################
      subroutine getq (iseq,ds,s0,h,q,valid,lfl)
c
c-----------------------------------------------------------------------
c
c ****** Obtain Q at the point given by the spherical coordinates in
c ****** vector S0 by tracing 5 field lines forwards and backwards
c ****** from S0 to the boundaries.
c
c ****** If Q was successfully obtained, return the Q value in
c ****** variable Q, and VALID=.T.; otherwise, return VALID=.F..
c
c ****** The length of the central field line is returned in LFL.
c
c ****** ISEQ is a sequence number that is used to construct the
c ****** file names when field line traces are being written to
c ****** output files.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use field
      use debug
      use diags
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iseq
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0
      real(r_typ) :: h
      real(r_typ) :: q
      logical :: valid
      real(r_typ) :: lfl
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Set the limit for how to evaluate Q derivatives near the
c ****** poles.
c
      real(r_typ), parameter :: st_pole_limit_max=5.e-3_r_typ
c
c ****** When sin(theta) of the central field line is smaller than
c ****** ST_POLE_LIMIT_MAX, Cartesian basis vectors are used to
c ****** compute derivatives of the field line mapping; otherwise,
c ****** spherical basis vectors are used.
c
c-----------------------------------------------------------------------
c
c ****** Flag to write warning messages.
c
      logical, parameter :: write_warning_messages=.true.
c
c-----------------------------------------------------------------------
c
      logical :: ttb
      type(flparam) :: ds_local
      real(r_typ), dimension(3) :: x0,x00
      real(r_typ), dimension(3) :: s00
      real(r_typ), dimension(3) :: bs0,bc0
      real(r_typ), dimension(3) :: bs1_f,bs1_b
      real(r_typ) :: s
      real(r_typ) :: bmag,br0_f,br0_b
      real(r_typ), dimension(3) :: bhat,bhat_abs
      real(r_typ) :: r_f,r_b,t_f,t_b
      real(r_typ) :: st_f,st_b
      integer, dimension(1) :: index_min
      real(r_typ), dimension(3) :: e,e1,e2
      real(r_typ), dimension(3) :: s1_f,s1_b
      real(r_typ), dimension(3) :: s1_1p_f,s1_1m_f
      real(r_typ), dimension(3) :: s1_1p_b,s1_1m_b
      real(r_typ), dimension(3) :: s1_2p_f,s1_2m_f
      real(r_typ), dimension(3) :: s1_2p_b,s1_2m_b
      real(r_typ) :: a_f,b_f,c_f,d_f
      real(r_typ) :: a_b,b_b,c_b,d_b
      real(r_typ) :: nsq
      type(csvec) :: xcs
      type(inout) :: outside
      logical :: outside_1p,outside_1m
      logical :: outside_2p,outside_2m
      real(r_typ) :: dx_1p,dx_1m,dx_2p,dx_2m
      real(r_typ) :: x_1p_f,x_1m_f,y_1p_f,y_1m_f
      real(r_typ) :: x_2p_f,x_2m_f,y_2p_f,y_2m_f
      real(r_typ) :: x_1p_b,x_1m_b,y_1p_b,y_1m_b
      real(r_typ) :: x_2p_b,x_2m_b,y_2p_b,y_2m_b
      character(6) :: ch_seq
      logical :: save_trace_points
      logical :: forward_ttb,backward_ttb
c
c ****** Field line trace storage buffer.
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: modulo_twopi
      logical, external :: outside_domain
c
c-----------------------------------------------------------------------
c
c ****** Set the flag to save field line traces.
c
      if (debug_level.ge.2) then
        save_trace_points=.true.
        call allocate_trajectory_buffer (xt)
        write (ch_seq,'(i6.6)') iseq
      else
        save_trace_points=.false.
      end if
c
c ****** Initialize Q and the validity flag.

      q=0.
      valid=.false.
      lfl=0.
c
c ****** Initilaize the local DS from that supplied in the
c ****** argument list.
c
      ds_local=ds
c
c ****** Set the flag to interpret the tracing direction as the
c ****** direction along the magnetic field line.
c
      ds_local%direction_is_along_b=.true.
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### COMMENT from GETQ:'
        write (*,*) '### Diagnostics for computation of Q:'
        write (*,*) 's0=',s0
      end if
c
c ****** Check that the requested position is not outside the domain.
c ****** If it is, return without calculating a valid Q.
c
      xcs%s=s0
      call sph_to_cart (xcs)
c
      if (outside_domain(b,xcs,outside)) return
c
c-----------------------------------------------------------------------
c ****** Trace the central field line.
c-----------------------------------------------------------------------
c
      s00=s0
c
      ds_local%direction=1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb)
      end if
c
      if (debug_level.ge.2) then
c$omp critical
        call write_trace ('fl_00_f_'//ch_seq//'.dat',xt)
c$omp end critical
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, forward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_f
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_f
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
c
      forward_ttb=ttb
c
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//
     &                ' the domain boundaries'
          write (*,*) '### during the forward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      else
        lfl=lfl+s
      end if
c
      ds_local%direction=-1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb)
      end if
c
      if (debug_level.ge.2) then
c$omp critical
        call write_trace ('fl_00_b_'//ch_seq//'.dat',xt)
c$omp end critical
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, backward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_b
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_b
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
c
      backward_ttb=ttb
c
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//
     &                ' the domain boundaries'
          write (*,*) '### during the backward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      else
        lfl=lfl+s
      end if
c
c ****** Save the central field line endpoints.
c
      r_f=s1_f(1)
      r_b=s1_b(1)
      t_f=s1_f(2)
      t_b=s1_b(2)
      st_f=sin(t_f)
      st_b=sin(t_b)
c
c ****** Get the radial component of the magnetic field at the
c ****** central field line endpoints.
c
      br0_f=bs1_f(1)
      br0_b=bs1_b(1)
c
c ****** If the field line did not reach both boundaries,
c ****** set the field line length to 0.
c
      if (.not.(forward_ttb.and.backward_ttb)) then
        lfl=0.
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line:'
        write (*,*) 'r_f=',r_f
        write (*,*) 'r_b=',r_b
        write (*,*) 't_f=',t_f
        write (*,*) 't_b=',t_b
        write (*,*) 'br0_f=',br0_f
        write (*,*) 'br0_b=',br0_b
        write (*,*) 'length=',lfl
      end if
c
c ****** Generate the basis vectors that define the plane of
c ****** the Q computation at the starting location
c ****** (in Cartesian coordinates).
c
c ****** Initialize the position vector at the starting point
c ****** in Cartesian coordinates in X0.
c
      call s2c (s0,x0)
c
c ****** Get the Cartesian magnetic field vector at S0.
c
      call sv_to_cv (s0,bs0,bc0)
c
      bmag=sqrt(bc0(1)**2+bc0(2)**2+bc0(3)**2)
c
      if (debug_level.ge.2) then
        write (*,*) 'bs0=',bs0
        write (*,*) 'bc0=',bc0
        write (*,*) 'bmag=',bmag
      end if
c
c ****** If we hit a null point (B=0), exit with an error.
c
      if (bmag.eq.0.) then
        if (debug_level.ge.2) then
          write (*,*)
          write (*,*) '### WARNING in GETQ:'
          write (*,*) 'B = 0 at the launch point.'
          write (*,*) 'Exiting with an error ...'
        end if
        return
      end if
c
      bhat=bc0/bmag
c
      if (debug_level.ge.2) then
        write (*,*) 'bhat=',bhat
      end if
c
c ****** Select the unit vector that is most perpendicular to B.
c
      bhat_abs=abs(bhat)
      index_min=minloc(bhat_abs)
c
      if (debug_level.ge.2) then
        write (*,*) 'index_min=',index_min
      end if
c
      e=0.
      e(index_min(1))=one
c
      if (debug_level.ge.2) then
        write (*,*) 'e=',e
      end if
c
c ****** The triplet E1, E2, and BHAT form an orthogonal basis.
c
      call normalized_cross_product (e,bhat,e1)
      call normalized_cross_product (e1,bhat,e2)
c
      if (debug_level.ge.2) then
        write (*,*) 'e1=',e1
        write (*,*) 'e2=',e2
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0+H*E1.
c-----------------------------------------------------------------------
c
      x00=x0+h*e1
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_1p=.true.
        dx_1p=0.
        s1_1p_f=s1_f
        s1_1p_b=s1_b
      else
        outside_1p=.false.
        dx_1p=h
      end if
c
      if (.not.outside_1p) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_p0_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_p0_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0-H*E1.
c-----------------------------------------------------------------------
c
      x00=x0-h*e1
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_1m=.true.
        dx_1m=0.
        s1_1m_f=s1_f
        s1_1m_b=s1_b
      else
        outside_1m=.false.
        dx_1m=h
      end if
c
c ****** If both the plus and minus perturbed launch points
c ****** are outside the domain, return without calculating a
c ****** valid Q.  Write a warning if this happens.
c
      if (outside_1p.and.outside_1m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e1 and x0-h*e1 launch points'//
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e1: ',s00
          x00=x0-h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e1: ',s00
        end if
        return
      end if
c
      if (.not.outside_1m) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_m0_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_m0_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0+H*E2.
c-----------------------------------------------------------------------
c
      x00=x0+h*e2
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_2p=.true.
        dx_2p=0.
        s1_2p_f=s1_f
        s1_2p_b=s1_b
      else
        outside_2p=.false.
        dx_2p=h
      end if
c
      if (.not.outside_2p) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0p_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0p_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
c
      end if
c
c-----------------------------------------------------------------------
c ****** Trace the field line at X0-H*E2.
c-----------------------------------------------------------------------
c
      x00=x0-h*e2
c
      call c2s (x00,s00)
c
c ****** If the launch point is outside the domain, use the central
c ****** field line to take a one-sided derivative.
c
      xcs%s=s00
      xcs%c=x00
c
      if (outside_domain(b,xcs,outside)) then
        outside_2m=.true.
        dx_2m=0.
        s1_2m_f=s1_f
        s1_2m_b=s1_b
      else
        outside_2m=.false.
        dx_2m=h
      end if
c
c ****** If both the plus and minus perturbed launch points
c ****** are outside the domain, return without calculating a
c ****** valid Q.  This should never happen.
c
      if (outside_2p.and.outside_2m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e2 and x0-h*e2 launch points'//
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e2: ',s00
          x00=x0-h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e2: ',s00
        end if
        return
      end if
c
      if (.not.outside_2m) then
c
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0m_f_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
c
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb)
        end if
c
        if (debug_level.ge.3) then
c$omp critical
          call write_trace ('fl_0m_b_'//ch_seq//'.dat',xt)
c$omp end critical
        end if
c
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
c
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
c
      end if
c
c ****** Get the value of the Jacobian matrix coefficients.
c
      if (st_f.lt.st_pole_limit_max) then
        x_1p_f=r_f*sin(s1_1p_f(2))*cos(s1_1p_f(3))
        x_1m_f=r_f*sin(s1_1m_f(2))*cos(s1_1m_f(3))
        y_1p_f=r_f*sin(s1_1p_f(2))*sin(s1_1p_f(3))
        y_1m_f=r_f*sin(s1_1m_f(2))*sin(s1_1m_f(3))
        x_2p_f=r_f*sin(s1_2p_f(2))*cos(s1_2p_f(3))
        x_2m_f=r_f*sin(s1_2m_f(2))*cos(s1_2m_f(3))
        y_2p_f=r_f*sin(s1_2p_f(2))*sin(s1_2p_f(3))
        y_2m_f=r_f*sin(s1_2m_f(2))*sin(s1_2m_f(3))
        a_f=(x_1p_f-x_1m_f)/(dx_1p+dx_1m)
        b_f=(x_2p_f-x_2m_f)/(dx_2p+dx_2m)
        c_f=(y_1p_f-y_1m_f)/(dx_1p+dx_1m)
        d_f=(y_2p_f-y_2m_f)/(dx_2p+dx_2m)
      else
        a_f=modulo_twopi(s1_1p_f(3)-s1_1m_f(3))*r_f*st_f/(dx_1p+dx_1m)
        b_f=modulo_twopi(s1_2p_f(3)-s1_2m_f(3))*r_f*st_f/(dx_2p+dx_2m)
        c_f=(s1_1p_f(2)-s1_1m_f(2))*r_f/(dx_1p+dx_1m)
        d_f=(s1_2p_f(2)-s1_2m_f(2))*r_f/(dx_2p+dx_2m)
      end if
c
      if (st_b.lt.st_pole_limit_max) then
        x_1p_b=r_b*sin(s1_1p_b(2))*cos(s1_1p_b(3))
        x_1m_b=r_b*sin(s1_1m_b(2))*cos(s1_1m_b(3))
        y_1p_b=r_b*sin(s1_1p_b(2))*sin(s1_1p_b(3))
        y_1m_b=r_b*sin(s1_1m_b(2))*sin(s1_1m_b(3))
        x_2p_b=r_b*sin(s1_2p_b(2))*cos(s1_2p_b(3))
        x_2m_b=r_b*sin(s1_2m_b(2))*cos(s1_2m_b(3))
        y_2p_b=r_b*sin(s1_2p_b(2))*sin(s1_2p_b(3))
        y_2m_b=r_b*sin(s1_2m_b(2))*sin(s1_2m_b(3))
        a_b=(x_1p_b-x_1m_b)/(dx_1p+dx_1m)
        b_b=(x_2p_b-x_2m_b)/(dx_2p+dx_2m)
        c_b=(y_1p_b-y_1m_b)/(dx_1p+dx_1m)
        d_b=(y_2p_b-y_2m_b)/(dx_2p+dx_2m)
      else
        a_b=modulo_twopi(s1_1p_b(3)-s1_1m_b(3))*r_b*st_b/(dx_1p+dx_1m)
        b_b=modulo_twopi(s1_2p_b(3)-s1_2m_b(3))*r_b*st_b/(dx_2p+dx_2m)
        c_b=(s1_1p_b(2)-s1_1m_b(2))*r_b/(dx_1p+dx_1m)
        d_b=(s1_2p_b(2)-s1_2m_b(2))*r_b/(dx_2p+dx_2m)
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Jacobian matrix coefficients:'
        write (*,*) 'a_f=',a_f
        write (*,*) 'b_f=',b_f
        write (*,*) 'c_f=',c_f
        write (*,*) 'd_f=',d_f
        write (*,*) 'a_b=',a_b
        write (*,*) 'b_b=',b_b
        write (*,*) 'c_b=',c_b
        write (*,*) 'd_b=',d_b
      end if
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Div B conservation check:'
        write (*,*) 'abs(Det(D_f))*abs(br0_f)/bmag = ',
     &              abs(a_f*d_f-b_f*c_f)*abs(br0_f)/bmag
        write (*,*) 'abs(Det(D_b))*abs(br0_b)/bmag = ',
     &              abs(a_b*d_b-b_b*c_b)*abs(br0_b)/bmag
      end if
c
c ****** Get the value of Q.
c
      nsq= (a_f*d_b-b_f*c_b)**2
     &    +(a_b*b_f-a_f*b_b)**2
     &    +(d_b*c_f-d_f*c_b)**2
     &    +(a_b*d_f-b_b*c_f)**2
c
      q=nsq*abs(br0_f*br0_b)/bmag**2
c
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Final Q value:'
        write (*,*) 'q=',q
      end if
c
      valid=.true.
c
c ****** Deallocate the storage for the field line trace buffer
c ****** if it was used.
c
      if (save_trace_points) then
        call deallocate_trajectory_buffer (xt)
      end if
c
      return
      end
c#######################################################################
      subroutine normalized_cross_product (a,b,c)
c
c-----------------------------------------------------------------------
c
c ****** Return the unit vector C = (A x B)/|A x B|.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: a,b,c
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: cnorm
c
c-----------------------------------------------------------------------
c
c ****** Set C to the cross-product of A and B.
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
c
c ****** Normalize C to unit length.
c
      cnorm=sqrt(c(1)**2+c(2)**2+c(3)**2)
c
      if (cnorm.ne.0.) then
        c=c/cnorm
      end if
c
      return
      end
c#######################################################################
      subroutine write_trace (fname,xt)
c
c-----------------------------------------------------------------------
c
c ****** Write the field line trace in structure XT to the
c ****** text file named FNAME.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      character, parameter :: TAB=achar(9)
c
c-----------------------------------------------------------------------
c
      integer :: i,ierr
c
c-----------------------------------------------------------------------
c
      call ffopen (1,fname,'rw',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_TRACE:'
        write (*,*) '### Could not create a field line output'//
     &              ' file.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
c ****** Write the output coordinates.
c
      write (1,'(5a)') 'r',TAB,'t',TAB,'p'
      do i=1,xt%npts
        write (1,'(3(1pe23.16,a))') xt%x(1)%f(i),TAB,
     &                              xt%x(2)%f(i),TAB,
     &                              xt%x(3)%f(i)
      enddo
c
      close (1)
c
      return
      end
c#######################################################################
      subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,
     &                    traced_to_r_boundary,xt)
c
c-----------------------------------------------------------------------
c
c ****** Trace a magnetic field line.
c
c-----------------------------------------------------------------------
c
c ****** The 3D magnetic field is specified by structure B.
c ****** The structure DS has the field line integration parameters,
c ****** and S0 contains the spherical coordinates of the launch
c ****** point.
c
c ****** TRACED_TO_R_BOUNDARY=.T. is returned if the field line
c ****** was traced all the way to a radial domain boundary, in
c ****** which case S1 contains the spherical coordinates of the
c ****** final location, and S has the traced field line length.
c
c ****** Otherwise, TRACED_TO_R_BOUNDARY=.F. is returned, and
c ****** S and S1 do not necessarily have valid values.
c ****** This also occurs if |B|=0 is encountered during the trace.
c
c ****** The magnetic field vectors in spherical coordinates
c ****** at the starting and ending footpoints, respectively,
c ****** are returned in BS0 and BS1.
c
c ****** If the field line trace is needed, pass in the optional
c ****** field line trace buffer XT.  The buffer XT needs to be
c ****** allocated prior to being passed in to this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use step_size_stats
      use debug
      use integrate_fl
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0,s1
      real(r_typ), dimension(3) :: bs0,bs1
      real(r_typ) :: s
      logical :: traced_to_r_boundary
      type(traj), optional :: xt
c
      intent(in) :: b,ds,s0
      intent(out) :: s1,bs0,bs1,s,traced_to_r_boundary
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      logical :: store_trace
      real(r_typ) :: ds0,dss,dsss,frac,dsmult_corrector
      real(r_typ) :: sf=1._r_typ,sf1,sf2
      logical :: done_tracing,first,nullb
      integer :: idir0,n,ntry,max_n,max_ntry
      type (csvec) :: x,xp,xo,bv,bhat1,bhat2
      type (inout) :: outside
      integer :: ierr
      real(r_typ) :: arc_length_remaining,ds_b
      type(flparam) :: current_ds
      logical :: tfc
c
      integer :: local_stat_n
      real(r_typ) :: local_stat_ds_sum,
     &               local_stat_ds_min,
     &               local_stat_ds_max
c
c-----------------------------------------------------------------------
c
      logical, external :: outside_domain
c
c-----------------------------------------------------------------------
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### Starting a new trace:'
        write (*,*) 'S0 = ',S0
      end if
c
c ****** Initializations.
c
      store_trace=present(xt)
c
      current_ds=ds
      first=.true.
      tfc=.false.
      s1=s0
      max_n=floor(ds%lmax/ds%min)
      max_ntry=ds%short_fl_max_tries
c
      if (gather_stats.gt.0) then
        local_stat_n=0
        local_stat_ds_sum=0.
        local_stat_ds_min=huge(one)
        local_stat_ds_max=0.
      end if
c
c ****** If storage for the field line trace was requested,
c ****** initialize the field line trace buffer.
c
      if (store_trace) xt%npts=0
c
      do ntry=1,max_ntry
c
c-----------------------------------------------------------------------
c ****** Trace a field line starting at launch point S0.
c-----------------------------------------------------------------------
c
        done_tracing=.false.
c
c ****** If things go wrong, return TRACED_TO_R_BOUNDARY=.F..
c
        traced_to_r_boundary=.false.
c
        x%s=s0
        call sph_to_cart (x)
        s=0.
        if (store_trace) call add_trajectory_point (xt,x%s)
c
c ****** Set the starting step size.
c
        ds0=current_ds%min
c
        if (debug_level.ge.4) then
          write (*,*)
          write (*,*) 'Start of trace:'
          write (*,*) 'X = ',x%s
        end if
c
c ****** If the initial point is outside the domain, stop tracing.
c
        if (outside_domain(b,x,outside)) tfc=.true.

        do n=1,max_n
c
c ****** If file line is done tracing (including restarts), exit loop.
c
          if(tfc) exit
c
c-----------------------------------------------------------------------
c ****** Trace to the next step.
c-----------------------------------------------------------------------
c
c ****** Set the tracing direction.
c
          call getb (b,x,bv)
c
c ****** On the first time in, set the tracing direction in
c ****** variable IDIR0 based on that specified by DS%DIRECTION
c ****** and DS%DIRECTION_IS_ALONG_B.  Store the magnetic field at
c ****** the initial location (in spherical coordinates) in BS0.
c
c ****** When DS%DIRECTION_IS_ALONG_B=.F., DS%DIRECTION=1 traces
c ****** in the direction of increasing r, and DS%DIRECTION=-1
c ****** traces in the opposite direction.
c
c ****** When DS%DIRECTION_IS_ALONG_B=.T., DS%DIRECTION=1 traces
c ****** in the direction of the magnetic field vector, and
c ****** DS%DIRECTION=-1 traces in the opposite direction.
c
          if (first) then
            first=.false.
            bs0=bv%s
            if (ds%direction.gt.0) then
              idir0=1
            else
              idir0=-1
            end if
            if (.not.ds%direction_is_along_b) then
              if (bv%s(1).lt.0.) idir0=-idir0
            end if
          end if
c
          call normalize_v (bv,nullb)
          if (nullb) then
            write (*,*)
            write (*,*) '### WARNING from TRACEFL:'
            write (*,*) '### The trace encountered a null pnt (B = 0).'
            write (*,*) '### This occurred at the start of the trace.'
            write (*,*) '### Abandoning the trace ...'
            write (*,*) 'Location (r,t,p) = ',x%s
            tfc=.true.
            exit
          end if
c
          if (n.gt.1) then
            ds_b=abs(dsss)
            bhat1=bhat2
          end if
          bhat2=bv
c
c ****** If a variable step size is being used, set the step
c ****** size for the next step.
c
c ****** If this is a repeat trace (of a short field line,
c ****** NTRY.gt.1), then leave the step size at the minimum
c ****** value until the number of points exceeds the minimum
c ****** allowed number.
c
          if (ds%variable.and.n.gt.1) then
            if (.not.(ntry.gt.1.and.n.le.ds%short_fl_min_points)) then
              call get_ds (b,x,bhat1,bhat2,ds_b,current_ds,ds0)
            end if
          end if
c
c ****** Check to see if this trace segment ends the trace
c ****** (i.e., exceeds the arc length specified, DS%LMAX).
c
          arc_length_remaining=ds%lmax-s
c
          if (ds0.ge.arc_length_remaining) then
            dss=arc_length_remaining
            done_tracing=.true.
          else
            dss=ds0
          end if
c
          if (dss.le.0.) then
            tfc=.true.
            exit
          end if
c
c ****** Gather step size statistics.
c
          if (gather_stats.gt.0) then
            local_stat_n=local_stat_n+1
            local_stat_ds_sum=local_stat_ds_sum+abs(ds0)
            local_stat_ds_min=min(local_stat_ds_min,abs(ds0))
            local_stat_ds_max=max(local_stat_ds_max,abs(ds0))
          end if
c
c-----------------------------------------------------------------------
c ****** Predictor.
c-----------------------------------------------------------------------
c
c ****** Advance for a half-step to achieve second-order accuracy.
c
          dsmult_corrector=one
          xo=x
          xp=x
c
          if (debug_level.ge.4) then
            write (*,*)
            write (*,*) 'Predictor:'
            write (*,*) 'N = ',n
            write (*,*) 'BV = ',bv%s
          end if
c
          dsss=half*idir0*dss
          call advance (xp,bv,dsss)
c
          if (debug_level.ge.4) then
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'XP = ',xp%s
          end if
c
c ****** Check if the field line has exited the domain.
c
          if (outside_domain(b,xp,outside)) then
c
            if (outside%t.or.outside%p) then
c
c ****** The field line has crossed the theta or phi boundaries;
c ****** stop tracing.
c
              tfc=.true.
              exit
c
            else if (outside%r) then
c
c ****** The field line has crossed the r boundary.
c
              if (debug_level.ge.4) then
                write (*,*) '### Predictor: Outside r domain:'
              end if
c
c ****** Get the clip fraction, FRAC, that clips the segment
c ****** to the r boundary.
c
              call get_r_clip_fraction (b,xo,xp,outside%r0,frac,ierr)
c
              if (debug_level.ge.4) then
                write (*,*) 'After GET_R_CLIP_FRACTION (predictor):'
                write (*,*) 'IERR = ',ierr
                write (*,*) 'FRAC = ',frac
              end if
c
c ****** If there is an error in getting the clip fraction,
c ****** stop tracing.  This should never happen: write
c ****** detailed debugging information.
c
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Predictor, 1st step:'
                write (*,*) '### Could not get the clip fraction.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                      ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'XP = ',xp
                tfc=.true.
                exit
              end if
c
              if (frac.eq.0..and.n.eq.1) then
c
c ****** The initial point is exactly on the boundary, and is
c ****** being traced out of the boundary.  We are done.
c
                traced_to_r_boundary=.true.
                done_tracing=.true.
                tfc=.true.
                exit
              end if
c
              if (frac.le.ds%predictor_min_clip_fraction) then
c
c ****** The starting point is close to the boundary.  Clip the
c ****** final point to the radial boundary and stop tracing.
c
                call clip_to_r_boundary (b,xo,xp,outside%r0,frac)
c
                if (debug_level.ge.4) then
                  write (*,*) 'After CLIP_TO_R_BOUNDARY (predictor):'
                  write (*,*) 'XP = ',xp%s
                end if
c
                x=xp
                dsss=frac*dsss
                traced_to_r_boundary=.true.
                done_tracing=.true.
                if (do_integral_along_fl) then
                   call getsf (xo,sf1)
                   call getsf (xp,sf2)
                   sf=half*(sf1+sf2)
                endif
                s=s+abs(dsss)*sf
                exit
c
              end if
c
c ****** If the starting point is not close enough to the radial
c ****** boundary, predict again half way to the radial boundary
c ****** to achieve second-order accuracy in the corrector.
c
              dsss=half*frac*dsss
              xp=x
              call advance (xp,bv,dsss)
c
              if (debug_level.ge.4) then
                write (*,*) 'After 2nd predictor:'
                write (*,*) 'XP = ',xp%s
              end if
c
c ****** Check if the predicted point has exited the domain.
c ****** This should never happen: write detailed debugging
c ****** information.
c
              if (outside_domain(b,xp,outside)) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Predictor, 2nd step:'
                write (*,*) '### Point is outside the boundary.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                  ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'XP = ',xp
                tfc=.true.
                exit
              end if
c
c ****** Reduce the step size for the corrector.  This minimizes the
c ****** chance of "going through the radial boundary and back into
c ****** the domain again", which can happen when the step size
c ****** is too big.
c
              dsmult_corrector=half*frac
c
            end if
c
          end if
c
c-----------------------------------------------------------------------
c ****** Corrector.
c-----------------------------------------------------------------------
c
          call getb (b,xp,bv)
          call normalize_v (bv,nullb)
          if (nullb) then
            write (*,*)
            write (*,*) '### WARNING from TRACEFL:'
            write (*,*) '### The trace encountered a null pnt (B = 0).'
            write (*,*) '### This occurred during the corrector.'
            write (*,*) '### Abandoning the trace ...'
            write (*,*) 'Location (r,t,p) = ',xp%s
            exit
          end if
c
          dsss=idir0*dsmult_corrector*dss
          call advance (x,bv,dsss)
c
          if (debug_level.ge.4) then
            write (*,*) 'After corrector advance:'
            write (*,*) 'BV = ',bv%s
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'X = ',x%s
          end if
c
          if (outside_domain(b,x,outside)) then
            if (outside%t.or.outside%p) then
              tfc=.true.
              exit
            else if (outside%r) then
c
              if (debug_level.ge.4) then
                write (*,*) '### Corrector: Outside r domain:'
              end if
c
              call get_r_clip_fraction (b,xo,x,outside%r0,frac,ierr)
c
c ****** If there is an error in getting the clip fraction,
c ****** stop tracing.  This should never happen: write
c ****** detailed debugging information.
c
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Corrector step:'
                write (*,*) '### Could not get the clip fraction.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',
     &                      ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'X  = ',x
                tfc=.true.
                exit
              end if
c
              if (debug_level.ge.4) then
                write (*,*) 'After GET_R_CLIP_FRACTION (corrector):'
                write (*,*) 'IERR = ',ierr
                write (*,*) 'FRAC = ',frac
              end if
c
              done_tracing=.true.
              traced_to_r_boundary=.true.
c
              if (frac.eq.0.) then
c
c ****** This should only happen when the previous point was exactly
c ****** on the boundary, and is being traced out of the domain.
c ****** In this case, do not store this point: we are done.
c
                x=xo
c
                if (debug_level.ge.4) then
                  write (*,*) 'The trace went from the boundary to'//
     &                        ' the outside (corrector):'
                  write (*,*) 'X = ',x%s
                end if
c
                exit
c
              else
c
                call clip_to_r_boundary (b,xo,x,outside%r0,frac)
                dsss=frac*dsss
c
                if (debug_level.ge.4) then
                  write (*,*) 'After CLIP_TO_R_BOUNDARY (corrector):'
                  write (*,*) 'X = ',x%s
                end if
c
              end if
            end if
          end if
c
c ****** Increment the number of points and the arc length.
c
          if (do_integral_along_fl) then
            call getsf (xo,sf1)
            call getsf (x,sf2)
            sf=half*(sf1+sf2)
          endif
          s=s+abs(dsss)*sf
c
c ****** Add the current position to the field line buffer
c ****** if requested.
c
          if (store_trace) call add_trajectory_point (xt,x%s)
c
c ****** Break out of max_n loop if done.
c
          if (done_tracing) exit
c
        enddo !n->max_n
c
c ****** Finished tracing the field line.
c
c ****** Break out of outer loop if field line complete.
c
        if (tfc) exit
c
c ****** Check that the number of points in the field line trace
c ****** exceeds the minimum allowed, DS%SHORT_FL_MIN_POINTS.
c ****** If not, reduce the step size and retrace the field line,
c ****** up to a maximum of DS%SHORT_FL_MAX_TRIES times.
c
        if (n.lt.ds%short_fl_min_points) then
c
          if (debug_level.ge.4) then
            write (*,*) 'Short field line:'
            write (*,*) 'N = ',n
            write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
            write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
            write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
          end if
c
          current_ds%over_rc=current_ds%over_rc*
     &                       ds%short_fl_shrink_factor
          current_ds%min=current_ds%min*
     &                   ds%short_fl_shrink_factor
          current_ds%max=current_ds%max*
     &                   ds%short_fl_shrink_factor
c
        else
          exit !Break out of max_ntry loop if line had enough points.
        end if
c
      end do !ntry->max_ntry
c
      if((ntry.ge.max_ntry) .and. (n.le.ds%short_fl_min_points)) then
        write (*,*)
        write (*,*) '### WARNING from TRACEFL:'
        write (*,*) 'Short field line after max tries:'
        write (*,*) 'Number of points in field line = ',n
        write (*,*) 'Number of tries = ',ntry
        write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
        write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
        write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
        write (*,*) 'S0 = ',s0
        write (*,*) 'S1 = ',x%s
      end if
c
c ****** Update the step size statistics.
c
      if (gather_stats.gt.0) then
c$omp critical (omp_stat)
        stat_n=stat_n+local_stat_n
        stat_ds_sum=stat_ds_sum+local_stat_ds_sum
        stat_ds_min=min(stat_ds_min,local_stat_ds_min)
        stat_ds_max=max(stat_ds_max,local_stat_ds_max)
c$omp end critical (omp_stat)
      end if
c
c ****** Store the final location in S1.
c
      s1=x%s
c
c ****** If the field line was not traced to the radial boundary,
c ****** do not attempt to return the magnetic field vector
c ****** at the field line endpoint.
c
      if (.not.traced_to_r_boundary) then
        bs1=0.
        return
      end if
c
c ****** Store the magnetic field vector (in spherical coordinates)
c ****** at the field line endpoint in BS1.
c
      call getb (b,x,bv)
      bs1=bv%s
c
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### About to exit:'
        write (*,*) 'N = ',n
        write (*,*) 'S = ',s
        write (*,*) 'S1 = ',s1
        write (*,*) 'BS0 = ',bs0
        write (*,*) 'BS1 = ',bs1
      end if
c
      return
      end
c#######################################################################
      subroutine get_ds (b,x,v0,v1,ds_v,ds,deltas)
c
c-----------------------------------------------------------------------
c
c ****** Set the integration step size based on the radius of
c ****** curvature of the field line, as estimated from the unit
c ****** vectors V0 and V1, and the local mesh cell size of the
c ****** magnetic field B at X (if requested).
c
c ****** The vectors V0 and V1 are assumed to have been evaluated
c ****** along the field line a distance DS_V apart.
c
c ****** On input, DELTAS should have the present value of the
c ****** step size.  On return, DELTAS is overwritten by the new
c ****** step size.
c
c ****** It is assumed that DS_V and DELTAS are positive.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x
      type(csvec) :: v0,v1
      real(r_typ) :: ds_v
      type(flparam) :: ds
      real(r_typ) :: deltas
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dv,factor,ds_mesh
c
c-----------------------------------------------------------------------
c
c ****** If DELTAS is zero, set it to the minimum value.
c
      if (deltas.eq.0.) then
        deltas=ds%min
        return
      end if
c
c ****** First, set the step size based on the local radius
c ****** of curvature of the field line.
c
c ****** Compute the factor by which DELTAS must be multiplied to
c ****** achieve the specified ratio of step size to radius
c ****** of curvature, DS%OVER_RC.
c
      dv=sqrt( (v0%c(1)-v1%c(1))**2
     &        +(v0%c(2)-v1%c(2))**2
     &        +(v0%c(3)-v1%c(3))**2)
c
      if (dv.ne.0.) then
        factor=ds_v*ds%over_rc/(dv*deltas)
      else
        factor=huge(dv)
      end if
c
c ****** Limit the change in DELTAS to the maximum permitted
c ****** change per step.
c
      factor=min(factor,ds%max_increase_factor)
      factor=max(factor,ds%max_decrease_factor)
c
c ****** Set the new step size.
c
      deltas=factor*deltas
c
c ****** Next, if requested, limit DELTAS by the local mesh
c ****** cell size of B at X.
c
c ****** Note that this only needs to be done if DELTAS is bigger
c ****** than DS%MIN, since only then is there a possibility of
c ****** reducing DELTAS further.
c
      if (deltas.gt.ds%min.and.ds%limit_by_local_mesh) then
c
c ****** Get the local mesh cell size at X.
c
        call get_local_mesh_size (b,x%s,ds_mesh)
c
c ****** Set DELTAS so that it does not exceed DS%LOCAL_MESH_FACTOR
c ****** times the local mesh size.
c
        deltas=min(deltas,ds%local_mesh_factor*ds_mesh)
c
      end if
c
c ****** Limit DELTAS by the maximum and mimimum allowed values.
c
      deltas=max(deltas,ds%min)
      deltas=min(deltas,ds%max)
c
      return
      end
c#######################################################################
      subroutine get_local_mesh_size (b,s,ds)
c
c-----------------------------------------------------------------------
c
c ****** Get the local mesh size DS from the magnetic field in
c ****** structure B at the spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      real(r_typ), dimension(3) :: s
      real(r_typ) :: ds
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ar,at,ap,drv,dtv,dpv,stv
c
c-----------------------------------------------------------------------
c
c ****** Get the local size of the main mesh of B at the
c ****** specified point.
c
      call interp (b%nrs,b%rs,s(1),i,ip1,ar,b%rs_invtab)
      drv=(one-ar)*b%drs(i)+ar*b%drs(ip1)
c
      call interp (b%nts,b%ts,s(2),j,jp1,at,b%ts_invtab)
      dtv=(one-at)*b%dts(j)+at*b%dts(jp1)
      stv=(one-at)*b%sts(j)+at*b%sts(jp1)
c
      call interp (b%nps,b%ps,s(3),k,kp1,ap,b%ps_invtab)
      dpv=(one-ap)*b%dps(k)+ap*b%dps(kp1)
c
      ds=min(drv,s(1)*dtv,s(1)*stv*dpv)
c
      return
      end
c#######################################################################
      subroutine normalize_v (v,null)
c
c-----------------------------------------------------------------------
c
c ****** Normalize the vector V to return a unit vector along V.
c ****** If V has zero length, then set NULL=.T. and leave V
c ****** unchanged; otherwise, set NULL=.F..
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: v
      logical :: null
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: vmag
c
c-----------------------------------------------------------------------
c
c ****** Use the spherical representation to compute the norm.
c
      vmag=sqrt(v%s(1)**2+v%s(2)**2+v%s(3)**2)
c
      if (vmag.eq.0.) then
        null=.true.
      else
        null=.false.
        v%c=v%c/vmag
        v%s=v%s/vmag
      end if
c
      return
      end
c#######################################################################
      subroutine advance (x,v,ds)
c
c-----------------------------------------------------------------------
c
c ****** Advance the position vector X by the step DS using the
c ****** velocity vector V.
c
c ****** This routine updates both Cartesian and spherical
c ****** represenations of X in the dual represenations position
c ****** vector X.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x,v
      real(r_typ) :: ds
c
c-----------------------------------------------------------------------
c
c ****** Advance the Cartesian position.
c
      x%c=x%c+ds*v%c
c
c ****** Transform the new Cartesian position to spherical
c ****** coordinates.
c
      call cart_to_sph (x)
c
      return
      end
c#######################################################################
      subroutine get_r_clip_fraction (b,x0,x1,outside_r0,frac,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Get the fraction FRAC (between 0 and 1) that expresses
c ****** the normalized distance between X0 and X1
c ****** corresponding to the location of the radial boundary.
c
c ****** The boundary location is obtained from the 3D magnetic
c ****** field in structure B.
c
c ****** It is assumed that X0 and X1 lie on different sides
c ****** of the r boundary specified by flag OUTSIDE_R0.
c
c ****** FRAC can be used to clip the position to the
c ****** radial boundary.
c
c ****** For a normal return, IERR=0 is returned when it was possible
c ****** to estimate FRAC; in the case of an inconsistency,
c ****** IERR=1 is returned and FRAC is invalid.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rb,dr0,dr1,dssq,x0dotb
      real(r_typ) :: eps,term1,term2,ratio,disc
      integer :: ipm
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      if (outside_r0) then
c
c ****** The segment crossed the inner radial boundary.
c
        rb=b%lim0(1)
        ipm=-1
c
      else
c
c ****** The segment crossed the outer radial boundary.
c
        rb=b%lim1(1)
        ipm=1
c
      end if
c
c ****** Check that X0 and X1 are consistent with a radial
c ****** boundary crossing.
c
      dr0=rb-x0%s(1)
      dr1=rb-x1%s(1)
c
      if (dr0.eq.0..and.dr1.eq.0.) then
        ierr=1
        return
      end if
c
      if (dr0*dr1.gt.0.) then
        ierr=1
        return
      end if
c
c ****** Treat the special case when X0 is exactly on the boundary.
c
      if (dr0.eq.0.) then
        frac=0.
        return
      end if
c
      dssq= (x1%c(1)-x0%c(1))**2
     &     +(x1%c(2)-x0%c(2))**2
     &     +(x1%c(3)-x0%c(3))**2
c
      if (dssq.le.0.) then
        ierr=1
        return
      end if
c
      x0dotb= x0%c(1)*(x1%c(1)-x0%c(1))
     &       +x0%c(2)*(x1%c(2)-x0%c(2))
     &       +x0%c(3)*(x1%c(3)-x0%c(3))
c
c ****** Get the square root of the discriminant, expanding small
c ****** arguments for accuracy.
c
      eps=10*sqrt(spacing(one))
c
      term1=x0dotb**2
      term2=dssq*dr0*(two*x0%s(1)+dr0)
      if (term1+term2.lt.0.) then
        ierr=1
        return
      end if
      if (term1.ne.0.) then
        ratio=term2/term1
        if (abs(ratio).lt.eps) then
          disc=sqrt(term1)*(one+half*ratio)
        else
          disc=sqrt(term1+term2)
        end if
      else
        disc=sqrt(term1+term2)
      end if
c
      frac=(-x0dotb+ipm*disc)/dssq
c
      return
      end
c#######################################################################
      subroutine clip_to_r_boundary (b,x0,x1,outside_r0,frac)
c
c-----------------------------------------------------------------------
c
c ****** Reset the position X1 to the normalized distance FRAC
c ****** between X0 and X1.
c
c ****** The boundary location is obtained from the 3D magnetic
c ****** field in structure B.
c
c ****** The flag OUTSIDE_R0=.T. indicates that the lower
c ****** radial boundary was crossed in going from X0 to X1;
c ****** otherwise, the upper radial boundary was crossed.
c
c ****** In conjunction with GET_R_CLIP_FRACTION this routine
c ****** can be used to clip points that cross the radial boundaries
c ****** to the boundary.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rval
c
c-----------------------------------------------------------------------
c
c ****** Reset the Cartesian position in X1.
c
      x1%c=(one-frac)*x0%c+frac*x1%c
c
c ****** Transform the new Cartesian position to spherical
c ****** coordinates.
c
      call cart_to_sph (x1)
c
c ****** Set the radius to the appropriate radial boundary
c ****** value (exactly).  This is done to take care of roundoff.
c
      if (outside_r0) then
        rval=b%lim0(1)
      else
        rval=b%lim1(1)
      end if
c
      x1%s(1)=rval
c
c ****** Transform the new spherical position to Cartesian
c ****** coordinates.
c
      call sph_to_cart (x1)
c
      return
      end
c#######################################################################
      function outside_domain (b,x,outside)
c
c-----------------------------------------------------------------------
c
c ****** If the spherical position in X lies outside the limits
c ****** of the domain, return a function result of .T.;
c ****** otherwise, return .F..
c
c ****** The detailed in/out location for each coordinate
c ****** is set in structure OUTSIDE.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x
      type(inout) :: outside
      logical :: outside_domain
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: pv
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Get the value of phi in the main interval [0,2*pi].
c
      pv=fold(zero,twopi,x%s(3))
c
      outside%r0=x%s(1).lt.b%lim0(1)
      outside%t0=x%s(2).lt.b%lim0(2)
      outside%p0=pv.lt.b%lim0(3)
      outside%r1=x%s(1).gt.b%lim1(1)
      outside%t1=x%s(2).gt.b%lim1(2)
      outside%p1=pv.gt.b%lim1(3)
c
      outside%r=outside%r0.or.outside%r1
      outside%t=outside%t0.or.outside%t1
      outside%p=outside%p0.or.outside%p1
c
      outside%domain=outside%r.or.outside%t.or.outside%p
c
      outside_domain=outside%domain
c
      return
      end
c#######################################################################
      function fold (x0,x1,x)
c
c-----------------------------------------------------------------------
c
c ****** "Fold" X into the periodic interval [X0,X1].
c
c ****** On return, X is such that X0.le.X.lt.X1.
c
c-----------------------------------------------------------------------
c
c ****** It is assumed that X0 does not equal X1, as is physically
c ****** necessary.  If X0 and X1 are equal, the routine just
c ****** returns with FOLD=X.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xl
c
c-----------------------------------------------------------------------
c
      fold=x
c
      if (x0.eq.x1) return
c
      xl=x1-x0
c
      fold=mod(x-x0,xl)+x0
c
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
c
      return
      end
c#######################################################################
      subroutine cart_to_sph (x)
c
c-----------------------------------------------------------------------
c
c ****** Update the dual representation position vector X so that
c ****** the Cartesian position corresponds to the spherical position.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x
c
c-----------------------------------------------------------------------
c
      call c2s (x%c,x%s)
c
      return
      end
c#######################################################################
      subroutine sph_to_cart (x)
c
c-----------------------------------------------------------------------
c
c ****** Update the dual representation position vector X so that
c ****** the spherical position corresponds to the Cartesian position.
c
c-----------------------------------------------------------------------
c
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x
c
c-----------------------------------------------------------------------
c
      call s2c (x%s,x%c)
c
      return
      end
c#######################################################################
      subroutine c2s (x,s)
c
c-----------------------------------------------------------------------
c
c ****** Convert the vector X = (x,y,z) from Cartesian coordinates
c ****** to spherical coordinates S = (r,t,p).
c
c ****** This routine returns T and P in radians, in the
c ****** following range:
c
c          0. .le. t .le. pi
c          0. .le. p .lt. 2.*pi
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: x
      real(r_typ), dimension(3) :: s
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r,t,p
c
c-----------------------------------------------------------------------
c
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
c
      if (r.eq.0.) then
        t=0.
      else
        t=acos(x(3)/r)
      end if
c
      if (x(1).eq.0.) then
        if (x(2).ge.0.) then
          p= halfpi
        else
          p=-halfpi
        end if
      else
        p=atan2(x(2),x(1))
      end if
      if (p.lt.0.) p=p+twopi
c
      s(1)=r
      s(2)=t
      s(3)=p
c
      return
      end
c#######################################################################
      subroutine s2c (s,x)
c
c-----------------------------------------------------------------------
c
c ****** Convert the vector S = (r,t,p) from spherical coordinates
c ****** to Cartesian coordinates X = (x,y,z).
c
c ****** This routine assumes that T and P are in radians.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s
      real(r_typ), dimension(3) :: x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: rst
c
c-----------------------------------------------------------------------
c
      rst=s(1)*sin(s(2))
      x(1)=rst*cos(s(3))
      x(2)=rst*sin(s(3))
      x(3)=s(1)*cos(s(2))
c
      return
      end
c#######################################################################
      subroutine cv_to_sv (s,cv,sv)
c
c-----------------------------------------------------------------------
c
c ****** Convert the Cartesian vector CV to spherical vector SV
c ****** at spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s,cv,sv
      intent(in) :: s,cv
      intent(out) :: sv
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: st,ct,sp,cp
c
c-----------------------------------------------------------------------
c
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
c
      sv(1)= cv(1)*st*cp+cv(2)*st*sp+cv(3)*ct
      sv(2)= cv(1)*ct*cp+cv(2)*ct*sp-cv(3)*st
      sv(3)=-cv(1)*sp   +cv(2)*cp
c
      return
      end
c#######################################################################
      subroutine sv_to_cv (s,sv,cv)
c
c-----------------------------------------------------------------------
c
c ****** Convert the spherical vector SV to Cartesian vector CV
c ****** at spherical position S.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: s,sv,cv
      intent(in) :: s,sv
      intent(out) :: cv
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: st,ct,sp,cp
c
c-----------------------------------------------------------------------
c
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
c
      cv(1)= sv(1)*st*cp+sv(2)*ct*cp-sv(3)*sp
      cv(2)= sv(1)*st*sp+sv(2)*ct*sp+sv(3)*cp
      cv(3)= sv(1)*ct   -sv(2)*st
c
      return
      end
c#######################################################################
      subroutine getb (b,x,bv)
c
c-----------------------------------------------------------------------
c
c ****** Get the interpolated magnetic field vector BV at the
c ****** dual representation position vector X from the
c ****** magnetic field vector in structure B.
c
c ****** When USE_ANALYTIC_FUNCTION=.TRUE., get the magnetic
c ****** field components by calling a function rather
c ****** than using B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use evaluate_spline_3d_interface
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vec) :: b
      type(csvec) :: x,bv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: pv
      real(r_typ) :: br,bt,bp,st,ct,sp,cp
      real(r_typ), dimension(3) :: xs
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Get the value of phi in the main interval [0,2*pi].
c
      pv=fold(zero,twopi,x%s(3))
c
c ****** If we are using an analytic function to define the
c ****** magnetic field, call it.
c
      if (use_analytic_function) then
        xs(1)=x%s(1)
        xs(2)=x%s(2)
        xs(3)=pv
        call magnetic_field_function (0._r_typ,.true.,xs,bv%s)
        go to 100
      end if
c
c ****** Get Br, Bt, and Bp.
c
      if (b%cubic) then
c
c ****** Cubic spline interpolation.
c
        br=evaluate_spline_3d(b%spl%r,x%s(1),x%s(2),pv,
     &                        b%inv(1)%c(1),
     &                        b%inv(1)%c(2),
     &                        b%inv(1)%c(3))
c
        bt=evaluate_spline_3d(b%spl%t,x%s(1),x%s(2),pv,
     &                        b%inv(2)%c(1),
     &                        b%inv(2)%c(2),
     &                        b%inv(2)%c(3))
c
        bp=evaluate_spline_3d(b%spl%p,x%s(1),x%s(2),pv,
     &                        b%inv(3)%c(1),
     &                        b%inv(3)%c(2),
     &                        b%inv(3)%c(3))
c
      else
c
c ****** Linear interpolation.
c
        call interp_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),
     &                  b%r%scales(1)%f,
     &                  b%r%scales(2)%f,
     &                  b%r%scales(3)%f,
     &                  b%inv(1),
     &                  b%r%f,x%s(1),x%s(2),pv,br)
c
        call interp_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),
     &                  b%t%scales(1)%f,
     &                  b%t%scales(2)%f,
     &                  b%t%scales(3)%f,
     &                  b%inv(2),
     &                  b%t%f,x%s(1),x%s(2),pv,bt)
c
        call interp_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),
     &                  b%p%scales(1)%f,
     &                  b%p%scales(2)%f,
     &                  b%p%scales(3)%f,
     &                  b%inv(3),
     &                  b%p%f,x%s(1),x%s(2),pv,bp)
c
      end if
c
      bv%s(1)=br
      bv%s(2)=bt
      bv%s(3)=bp
c
  100 continue
c
c ****** Transform the spherical components of the magnetic field
c ****** to the Cartesian components.
c
      st=sin(x%s(2))
      ct=cos(x%s(2))
      sp=sin(pv)
      cp=cos(pv)
      bv%c(1)=bv%s(1)*st*cp+bv%s(2)*ct*cp-bv%s(3)*sp
      bv%c(2)=bv%s(1)*st*sp+bv%s(2)*ct*sp+bv%s(3)*cp
      bv%c(3)=bv%s(1)*ct   -bv%s(2)*st
c
      return
      end
c#######################################################################
      subroutine getsf (x,sf)
c
c-----------------------------------------------------------------------
c
c ****** Get the interpolated scalar field SF at the
c ****** dual representation position vector X from the
c ****** field SCALAR_FIELD.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use vars
      use integrate_fl
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(csvec) :: x
      real(r_typ) :: sf
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: pv
      real(r_typ), dimension(3) :: xs
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Get the value of phi in the main interval [0,2*pi].
c
      pv=fold(zero,twopi,x%s(3))
c
c ****** Linear interpolation.
c
        call interp_3d (scalar_field%dims(1),
     &                  scalar_field%dims(2),
     &                  scalar_field%dims(3),
     &                  scalar_field%scales(1)%f,
     &                  scalar_field%scales(2)%f,
     &                  scalar_field%scales(3)%f,
     &                  inv_sf,
     &                  scalar_field%f,x%s(1),x%s(2),pv,sf)

c
      return
      end
c#######################################################################
      subroutine interp_3d (nx,ny,nz,x,y,z,inv,f,xv,yv,zv,fv)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the 3D field FV at (XV,YV,ZV) from
c ****** array F(NX,NY,NZ), defined on the mesh X(NX) x Y(NY) x Z(NZ).
c ****** The structure INV holds the inverse interpolation tables.
c
c ****** Note that if the point (XV,YV,ZV) is outside the bounds of
c ****** the X x Y x Z mesh, FV=0. is returned.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny,nz
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ) :: xv,yv,zv,fv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ax,ay,az
c
c-----------------------------------------------------------------------
c
c ****** If the point is outside the data limits, return a
c ****** zero value.
c
      if (xv.lt.x(1).or.xv.gt.x(nx).or.
     &    yv.lt.y(1).or.yv.gt.y(ny).or.
     &    zv.lt.z(1).or.zv.gt.z(nz)) then
        fv=0.
        return
      end if
c
      call interp (nx,x,xv,i,ip1,ax,inv%c(1))
      call interp (ny,y,yv,j,jp1,ay,inv%c(2))
      call interp (nz,z,zv,k,kp1,az,inv%c(3))
c
      fv= (one-ax)*( (one-ay)*( (one-az)*f(i  ,j  ,k  )
     &                         +     az *f(i  ,j  ,kp1))
     &              +     ay *( (one-az)*f(i  ,jp1,k  )
     &                         +     az *f(i  ,jp1,kp1)))
     &   +     ax *( (one-ay)*( (one-az)*f(ip1,j  ,k  )
     &                         +     az *f(ip1,j  ,kp1))
     &              +     ay *( (one-az)*f(ip1,jp1,k  )
     &                         +     az *f(ip1,jp1,kp1)))
c
      return
      end
c#######################################################################
      subroutine allocate_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Allocate the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      xt%ndim=3
      xt%size=xt%initial_size
c
c ****** Allocate storage for the trajectory buffer.
c
      allocate (xt%x(xt%ndim))
c
      do i=1,xt%ndim
        allocate (xt%x(i)%f(xt%size))
      enddo
c
c ****** Initialize the current number of points in the
c ****** trajectory buffer.
c
      xt%npts=0
c
      return
      end
c#######################################################################
      subroutine deallocate_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Deallocate the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      if (.not.associated(xt%x)) return
c
      do i=1,xt%ndim
        if (associated(xt%x(i)%f)) then
          deallocate (xt%x(i)%f)
        end if
      enddo
c
      deallocate (xt%x)
c
      xt%npts=0
c
      return
      end
c#######################################################################
      subroutine add_trajectory_point (xt,x)
c
c-----------------------------------------------------------------------
c
c ****** Add the position vector X to the trajectory buffer XT.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
      real(r_typ), dimension(xt%ndim) :: x
c
c-----------------------------------------------------------------------
c
      integer :: i,n
c
c-----------------------------------------------------------------------
c
c ****** Increment the number of points in the trajectory.
c ****** If the buffer is full, expand it.
c
      n=xt%npts
      n=n+1
      if (n.gt.xt%size) call expand_trajectory_buffer (xt)
c
c ****** Add the point to the buffer.
c
      do i=1,xt%ndim
        xt%x(i)%f(n)=x(i)
      enddo
c
      xt%npts=n
c
      return
      end
c#######################################################################
      subroutine expand_trajectory_buffer (xt)
c
c-----------------------------------------------------------------------
c
c ****** Expand the trajectory buffer XT by doubling the number
c ****** of points in the buffer.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use debug
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(traj) :: xt
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:), pointer :: f
      integer :: i,n
c
c-----------------------------------------------------------------------
c
c ****** Double the current buffer size, and copy the current
c ****** contents to the expanded buffer.
c
      n=xt%size
      do i=1,xt%ndim
        allocate (f(2*n))
        f(1:n)=xt%x(i)%f(1:n)
        deallocate (xt%x(i)%f)
        xt%x(i)%f=>f
      enddo
      xt%size=2*n
      if (debug_level.ge.5) then
        write (*,*) 'Expanded a trajectory buffer to ',xt%size
      end if
c
      return
      end
c#######################################################################
      subroutine build_inverse_tables (s,inv)
c
c-----------------------------------------------------------------------
c
c ****** Build the inverse interpolation tables INV for the SDS
c ****** in structure S.
c
c ****** These arrays are used to to increase the efficiency
c ****** of interpolation lookups.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      type(vtab) :: inv
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Use a number of points for the inverse interpolation table
c ****** equal to the number in the original scale.
c
      do i=1,s%ndim
        inv%c(i)%n=s%dims(i)
        allocate (inv%c(i)%f(inv%c(i)%n))
        call getinv (s%scales(i)%f,s%dims(i),inv%c(i))
      enddo
c
      return
      end
c#######################################################################
      subroutine getinv (x,n,tab)
c
c-----------------------------------------------------------------------
c
c ****** Build an inverse interpolation table to increase the
c ****** efficiency of table look-up in a nonuniform mesh.
c
c ****** On input, the table X(N) is specified, together with the
c ****** number of points to use in the inverse interpolation
c ****** table, NU.
c
c ****** The output is a structure TAB with the inverse interpolation
c ****** table.  This structure has the following components:
c
c ******    N:  the number of points in the table;
c ******    D:  the inverse of the uniform table spacing;
c ******    F:  the inverse interpolation table.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
      use interp_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,k,ip1
      real(r_typ) :: dx,xv,alpha,en
c
c-----------------------------------------------------------------------
c
c ****** Check that the number of points is valid.
c
      if (tab%n.le.1) then
        write (*,*)
        write (*,*) '### ERROR in GETINV:'
        write (*,*) '### Invalid number of points specified'//
     &              ' for the inverse interpolation table.'
        write (*,*)
        write (*,*) 'Number of points = ',tab%n
        call exit (1)
      end if
c
c ****** Set the uniform interval to be used in the inverse
c ****** interpolation.
c
      dx=(x(n)-x(1))/(tab%n-one)
c
      if (dx.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in GETINV:'
        write (*,*) '### Invalid interval for the inverse'//
     &              ' interpolation table.'
        write (*,*)
        write (*,*) 'Interval = ',dx
        call exit (1)
      end if
c
      tab%d=one/dx
c
c ****** Build the inverse interpolation table.
c
      en=n
c
      do k=1,tab%n
        xv=x(1)+(k-one)*dx
        xv=max(xv,x(1))
        xv=min(xv,x(n))
        call interp (n,x,xv,i,ip1,alpha)
        tab%f(k)=i+alpha
        tab%f(k)=max(tab%f(k),one)
        tab%f(k)=min(tab%f(k),en)
      enddo
c
      return
      end
c#######################################################################
      subroutine interp (n,x,xv,i,ip1,alpha,tab)
c
c-----------------------------------------------------------------------
c
c ****** Find the interval I in table X(i), i=1,2,...,N, that encloses
c ****** the value XV, such that X(I).le.XV.le.X(I+1).
c ****** For the special case when N=1, XV must equal X(1) exactly.
c
c ****** This routine uses LOCATE_INTERVAL (from the SPLINE library)
c ****** to do the actual work.  If the interval is not found
c ****** LOCATE_INTERVAL terminates with an error.
c
c ****** This routine does not do the actual interpolation.  However,
c ****** the returned values of I, IP1 (which generally equals I+1),
c ****** and ALPHA can be used to get the interpolant.
c
c ****** The optional inverse interpolation table, TAB, can be
c ****** supplied to improve the efficiency of the search.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i
      integer :: ip1
      real(r_typ) :: alpha
      type(itab), optional :: tab
      intent(in) :: n,x,xv,tab
      intent(out) :: i,ip1,alpha
c
c-----------------------------------------------------------------------
c
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
c
      if (n.eq.1) then
        ip1=1
        alpha=0.
      else
        ip1=i+1
        if (x(i).eq.x(i+1)) then
          alpha=0.
        else
          alpha=(xv-x(i))/(x(i+1)-x(i))
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine compute_spline_1d (nx,x,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 1D function
c ****** defined by array F(NX) with scale X(NX).
c
c ****** The spline coefficients are returned in structure S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(nx) :: f
      type(spl1d) :: s
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
c
      allocate (s%x(nx))
      allocate (s%f(nx))
      allocate (s%fxx(nx))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%f=f
c
      call ezspline (nx,x,s%f,s%fxx)
c
      return
      end
c#######################################################################
      subroutine compute_spline_2d (nx,ny,x,y,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 2D function
c ****** defined by array F(NX,NY), with scales X(NX) and Y(NY).
c
c ****** The spline coefficients are returned in structure S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      type(spl2d) :: s
c
c-----------------------------------------------------------------------
c
      integer :: i,j
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
      s%ny=ny
c
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%f(nx,ny))
      allocate (s%fxx(nx,ny))
      allocate (s%fyy(nx,ny))
      allocate (s%fxxyy(nx,ny))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%y=y
      s%f=f
c
      do j=1,ny
        gx(:)=s%f(:,j)
        call ezspline (nx,x,gx,gppx)
        s%fxx(:,j)=gppx(:)
      enddo
c
      do i=1,nx
        gy(:)=s%f(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fyy(i,:)=gppy(:)
      enddo
c
      do i=1,nx
        gy(:)=s%fxx(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fxxyy(i,:)=gppy(:)
      enddo
c
      return
      end
c#######################################################################
      subroutine compute_spline_3d (nx,ny,nz,x,y,z,f,s)
c
c-----------------------------------------------------------------------
c
c ****** Find the cubic spline coefficients for the 3D function
c ****** defined by array F(NX,NY,NZ), with scales X(NX), Y(NY),
c ****** and Z(NZ).
c
c ****** The spline coefficients are returned in structure S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny,nz
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      real(r_typ), dimension(nx,ny,nz) :: f
      type(spl3d) :: s
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
      real(r_typ), dimension(nz) :: gz,gppz
c
c-----------------------------------------------------------------------
c
c ****** Allocate storage for the spline coefficients.
c
      s%nx=nx
      s%ny=ny
      s%nz=nz
c
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%z(nz))
      allocate (s%f(nx,ny,nz))
      allocate (s%fxx(nx,ny,nz))
      allocate (s%fyy(nx,ny,nz))
      allocate (s%fzz(nx,ny,nz))
      allocate (s%fxxyy(nx,ny,nz))
      allocate (s%fxxzz(nx,ny,nz))
      allocate (s%fyyzz(nx,ny,nz))
      allocate (s%fxxyyzz(nx,ny,nz))
c
c ****** Evaluate the spline coefficients.
c
      s%x=x
      s%y=y
      s%z=z
      s%f=f
c
c$omp parallel default(shared)
c$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
c
c$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do j=1,ny
          gx(:)=s%f(:,j,k)
          call ezspline (nx,x,gx,gppx)
          s%fxx(:,j,k)=gppx(:)
        enddo
      enddo
c$omp end do
c
c$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do i=1,nx
          gy(:)=s%f(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fyy(i,:,k)=gppy(:)
        enddo
      enddo
c$omp end do
c
c$omp do collapse(2) schedule(dynamic)
      do j=1,ny
        do i=1,nx
          gz(:)=s%f(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fzz(i,j,:)=gppz(:)
        enddo
      enddo
c$omp end do
c$omp end parallel
c
c$omp parallel default(shared)
c$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
c$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do i=1,nx
          gy(:)=s%fxx(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fxxyy(i,:,k)=gppy(:)
        enddo
      enddo
c$omp end do
c$omp end parallel
c
c$omp parallel default(shared)
c$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
c$omp do collapse(2) schedule(dynamic)
      do j=1,ny
        do i=1,nx
          gz(:)=s%fxx(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxzz(i,j,:)=gppz(:)
          gz(:)=s%fyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fyyzz(i,j,:)=gppz(:)
          gz(:)=s%fxxyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxyyzz(i,j,:)=gppz(:)
        enddo
      enddo
c$omp end do
c$omp end parallel
c
      return
      end
c#######################################################################
      function evaluate_spline_1d (s,x,tab)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 1D spline in structure S at the
c ****** point X.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval that contains X.
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_1D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl1d) :: s
      real(r_typ) :: x
      type(itab), optional :: tab
      real(r_typ) :: evaluate_spline_1d
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
c
c-----------------------------------------------------------------------
c
c ****** Find the index of the cell enclosing point X.
c
      if (present(tab)) then
        i=locate_interval(s%nx,s%x,x,tab)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
c ****** Interpolate in x.
c
      evaluate_spline_1d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          s%f(i  ),
     &                          s%f(i+1),
     &                          s%fxx(i  ),
     &                          s%fxx(i+1),
     &                          x)
c
      return
      end
c#######################################################################
      function evaluate_spline_2d (s,x,y,tabx,taby)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 2D spline in structure S at the
c ****** point (X,Y).
c
c ****** The optional arguments TABX and TABY are inverse
c ****** interpolation tables that can be used to speed up the
c ****** search for the cell that contains (X,Y).
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_2D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl2d) :: s
      real(r_typ) :: x,y
      type(itab), optional :: tabx,taby
      real(r_typ) :: evaluate_spline_2d
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ii
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
c
c-----------------------------------------------------------------------
c
c ****** Flag to use the more efficient version of this routine.
c
      logical, parameter :: use_faster=.true.
c
c-----------------------------------------------------------------------
c
c ****** Find the indices of the cell enclosing (X,Y).
c
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
c
      if (use_faster) go to 100
c
c ****** Leff efficient version.  This version
c ****** is slower but slightly easier to understand.
c
c ****** Interpolate in y.
c
      do ii=0,1
        f1(ii)=speval(s%y(j  ),
     &                s%y(j+1),
     &                s%f(i+ii,j  ),
     &                s%f(i+ii,j+1),
     &                s%fyy(i+ii,j  ),
     &                s%fyy(i+ii,j+1),
     &                y)
        fxx1(ii)=speval(s%y(j  ),
     &                  s%y(j+1),
     &                  s%fxx(i+ii,j  ),
     &                  s%fxx(i+ii,j+1),
     &                  s%fxxyy(i+ii,j  ),
     &                  s%fxxyy(i+ii,j+1),
     &                  y)
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_2d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
c
  100 continue
c
c ****** More efficient version.  This version
c ****** is slighty faster than the one above.
c
c ****** Interpolate in y.
c
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
c
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,
     &                     s%f(i+ii,j  ),
     &                     s%f(i+ii,j+1),
     &                     s%fyy(i+ii,j  ),
     &                     s%fyy(i+ii,j+1))
        fxx1(ii)=speval_abcd(a,b,c,d,
     &                       s%fxx(i+ii,j  ),
     &                       s%fxx(i+ii,j+1),
     &                       s%fxxyy(i+ii,j  ),
     &                       s%fxxyy(i+ii,j+1))
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_2d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
      end
c#######################################################################
      function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the 3D spline in structure S at the
c ****** point (X,Y,Z).
c
c ****** The optional arguments TABX, TABY, and TABZ are inverse
c ****** interpolation tables that can be used to speed up the
c ****** search for the cell that contains (X,Y,Z).
c
c ****** The cubic spline coefficients for the spline S can
c ****** be obtained using routine COMPUTE_SPLINE_3D.
c
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(spl3d) :: s
      real(r_typ) :: x,y,z
      type(itab), optional :: tabx,taby,tabz
      real(r_typ) :: evaluate_spline_3d
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ii,jj
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ), dimension(0:1,0:1) :: f2,fxx2,fyy2,fxxyy2
      real(r_typ) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
c
c-----------------------------------------------------------------------
c
c ****** Flag to use the more efficient version of this routine.
c
      logical, parameter :: use_faster=.true.
c
c-----------------------------------------------------------------------
c
c ****** Find the indices of the cell enclosing (X,Y,Z).
c
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
c
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
c
      if (present(tabz)) then
        k=locate_interval(s%nz,s%z,z,tabz)
      else
        k=locate_interval(s%nz,s%z,z)
      end if
c
      if (use_faster) go to 100
c
c ****** Leff efficient version.  This version
c ****** is slower but slightly easier to understand.
c
c ****** Interpolate in z.
c
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval(s%z(k  ),
     &                     s%z(k+1),
     &                     s%f(i+ii,j+jj,k  ),
     &                     s%f(i+ii,j+jj,k+1),
     &                     s%fzz(i+ii,j+jj,k  ),
     &                     s%fzz(i+ii,j+jj,k+1),
     &                     z)
          fxx2(ii,jj)=speval(s%z(k  ),
     &                       s%z(k+1),
     &                       s%fxx(i+ii,j+jj,k  ),
     &                       s%fxx(i+ii,j+jj,k+1),
     &                       s%fxxzz(i+ii,j+jj,k  ),
     &                       s%fxxzz(i+ii,j+jj,k+1),
     &                       z)
          fyy2(ii,jj)=speval(s%z(k  ),
     &                       s%z(k+1),
     &                       s%fyy(i+ii,j+jj,k  ),
     &                       s%fyy(i+ii,j+jj,k+1),
     &                       s%fyyzz(i+ii,j+jj,k  ),
     &                       s%fyyzz(i+ii,j+jj,k+1),
     &                       z)
          fxxyy2(ii,jj)=speval(s%z(k  ),
     &                         s%z(k+1),
     &                         s%fxxyy(i+ii,j+jj,k  ),
     &                         s%fxxyy(i+ii,j+jj,k+1),
     &                         s%fxxyyzz(i+ii,j+jj,k  ),
     &                         s%fxxyyzz(i+ii,j+jj,k+1),
     &                         z)
        enddo
      enddo
c
c ****** Interpolate in y.
c
      do ii=0,1
        f1(ii)=speval(s%y(j  ),
     &                s%y(j+1),
     &                f2(ii,0),
     &                f2(ii,1),
     &                fyy2(ii,0),
     &                fyy2(ii,1),
     &                y)
        fxx1(ii)=speval(s%y(j  ),
     &                  s%y(j+1),
     &                  fxx2(ii,0),
     &                  fxx2(ii,1),
     &                  fxxyy2(ii,0),
     &                  fxxyy2(ii,1),
     &                  y)
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_3d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
c
  100 continue
c
c ****** More efficient version.  This version
c ****** is about 25% faster than the one above.
c
c ****** Interpolate in z.
c
      call speval_get_abcd (s%z(k),s%z(k+1),z,a,b,c,d)
c
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval_abcd(a,b,c,d,
     &                          s%f(i+ii,j+jj,k  ),
     &                          s%f(i+ii,j+jj,k+1),
     &                          s%fzz(i+ii,j+jj,k  ),
     &                          s%fzz(i+ii,j+jj,k+1))
          fxx2(ii,jj)=speval_abcd(a,b,c,d,
     &                            s%fxx(i+ii,j+jj,k  ),
     &                            s%fxx(i+ii,j+jj,k+1),
     &                            s%fxxzz(i+ii,j+jj,k  ),
     &                            s%fxxzz(i+ii,j+jj,k+1))
          fyy2(ii,jj)=speval_abcd(a,b,c,d,
     &                            s%fyy(i+ii,j+jj,k  ),
     &                            s%fyy(i+ii,j+jj,k+1),
     &                            s%fyyzz(i+ii,j+jj,k  ),
     &                            s%fyyzz(i+ii,j+jj,k+1))
          fxxyy2(ii,jj)=speval_abcd(a,b,c,d,
     &                              s%fxxyy(i+ii,j+jj,k  ),
     &                              s%fxxyy(i+ii,j+jj,k+1),
     &                              s%fxxyyzz(i+ii,j+jj,k  ),
     &                              s%fxxyyzz(i+ii,j+jj,k+1))
        enddo
      enddo
c
c ****** Interpolate in y.
c
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
c
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,
     &                     f2(ii,0),
     &                     f2(ii,1),
     &                     fyy2(ii,0),
     &                     fyy2(ii,1))
        fxx1(ii)=speval_abcd(a,b,c,d,
     &                       fxx2(ii,0),
     &                       fxx2(ii,1),
     &                       fxxyy2(ii,0),
     &                       fxxyy2(ii,1))
      enddo
c
c ****** Interpolate in x.
c
      evaluate_spline_3d=speval(s%x(i  ),
     &                          s%x(i+1),
     &                          f1(0),
     &                          f1(1),
     &                          fxx1(0),
     &                          fxx1(1),
     &                          x)
c
      return
      end
c#######################################################################
      subroutine spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** The boundary conditions at the two ends are specified
c ****** by IBC0 and IBC1, and the coefficient arrays C0 and C1.
c
c ****** These are defined at x=X(1) according to the value
c ****** of IBC0 and C0.  (The conditions at x=X(N) are
c ****** specified similarly by IBC1 and C1.)
c
c        IBC0 = 1:  Set the second derivative at one end of the
c                   cell to equal the second derivative at the
c                   other end of the cell. [f''(1)=f''(2)]
c
c        IBC0 = 2:  Set the second derivative to zero.
c                   This corresponds to a "natural spline".
c                   [f''(1)=0.]
c
c        IBC0 = 3:  Set the first derivative to zero.
c                   [f'(1)=0.]
c
c        IBC0 = 4:  Set the second derivative to C0(1).
c                   [f''(1)=C0(1)]
c
c        IBC0 = 5:  Set the first derivative to C0(1).
c                   [f'(1)=C0(1)]
c
c        IBC0 = 6:  Set a linear combination of the first and
c                   second derivatives, according to C0(1),
c                   C0(2), and C0(3).
c                   [C0(2)*f'(1)+C0(3)*f''(1)=C0(1)]
c
c        IBC0 = 7:  Set a linear combination of the second
c                   derivatives at the left and right ends
c                   of the first cell.
c                   [C0(1)*f''(1)+C0(2)*f''(2)=C0(3)]
c
c        IBC0 = 8:  Set a linear combination of the first
c                   derivatives at the left and right ends
c                   of the first cell.
c                   [C0(1)*f'(1)+C0(2)*f'(2)=C0(3)]
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
c
      intent(in) :: n,x,f,ibc0,c0,ibc1,c1
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 3 points.
c
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the boundary condition at X(1).
c-----------------------------------------------------------------------
c
c ****** Set the unused value to zero.
c
      c(1)=0.
c
      select case (ibc0)
      case (1)
c
c ****** Second derivatives equal at the left and right ends
c ****** of the first cell.
c
        a(1)=one
        b(1)=-one
        fpp(1)=0.
c
      case (2)
c
c ****** Second derivative is zero ("natural splines").
c
        a(1)=one
        b(1)=0.
        fpp(1)=0.
c
      case (3)
c
c ****** First derivative is zero.
c
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx
c
      case (4)
c
c ****** Second derivative is specified.
c
        a(1)=one
        b(1)=0.
        fpp(1)=c0(1)
c
      case (5)
c
c ****** First derivative is specified.
c
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx-c0(1)
c
      case (6)
c
c ****** A combination of the first and second derivatives
c ****** is specified.
c
        if (c0(2).eq.0..and.c0(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//
     &                ' at X(1).'
          write (*,*) '### Boundary condition type IBC0 = 6.'
          write (*,*) '### It is illegal for both C0(2) and C0(3)'//
     &                ' to be zero.'
          call exit (1)
        end if
c
        if (c0(2).ne.0.) then
          dx=x(2)-x(1)
          a(1)=dx*third-c0(3)/c0(2)
          b(1)=dx*sixth
          fpp(1)=(f(2)-f(1))/dx-c0(1)/c0(2)
        else
          a(1)=one
          b(1)=0.
          fpp(1)=c0(1)/c0(3)
        end if
c
      case (7)
c
c ****** A linear combination of the second derivatives at the left
c ****** and right ends of the first cell is specified.
c
        a(1)=c0(1)
        b(1)=c0(2)
        fpp(1)=c0(3)
c
      case (8)
c
c ****** A linear combination of the first derivatives at the left
c ****** and right ends of the first cell is specified.
c
        dx=x(2)-x(1)
        a(1)=dx*(c0(1)/3-c0(2)/6)
        b(1)=dx*(c0(1)/6-c0(2)/3)
        fpp(1)=(c0(1)+c0(2))*(f(2)-f(1))/dx-c0(3)
c
      case default
c
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(1).'
        write (*,*) '### IBC0 is invalid.'
        write (*,*) 'IBC0 = ',ibc0
        call exit (1)
c
      end select
c
c-----------------------------------------------------------------------
c ****** Set the boundary condition at X(N).
c-----------------------------------------------------------------------
c
c ****** Set the unused value to zero.
c
      b(n)=0.
c
      select case (ibc1)
      case (1)
c
c ****** Second derivatives equal at the left and right ends
c ****** of the last cell.
c
        a(n)=one
        c(n)=-one
        fpp(n)=0.
c
      case (2)
c
c ****** Second derivative is zero ("natural splines").
c
        a(n)=one
        c(n)=0.
        fpp(n)=0.
c
      case (3)
c
c ****** First derivative is zero.
c
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=-(f(n)-f(n-1))/dx
c
      case (4)
c
c ****** Second derivative is specified.
c
        a(n)=one
        c(n)=0.
        fpp(n)=c1(1)
c
      case (5)
c
c ****** First derivative is specified.
c
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=c1(1)-(f(n)-f(n-1))/dx
c
      case (6)
c
c ****** A combination of the first and second derivatives
c ****** is specified.
c
        if (c1(2).eq.0..and.c1(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//
     &                ' at X(N).'
          write (*,*) '### Boundary condition type IBC1 = 6.'
          write (*,*) '### It is illegal for both C1(2) and C1(3)'//
     &                ' to be zero.'
          call exit (1)
        end if
c
        if (c1(2).ne.0.) then
          dx=x(n)-x(n-1)
          a(n)=dx*third+c1(3)/c1(2)
          c(n)=dx*sixth
          fpp(n)=c1(1)/c1(2)-(f(n)-f(n-1))/dx
        else
          a(n)=one
          c(n)=0.
          fpp(n)=c1(1)/c1(3)
        end if
c
      case (7)
c
c ****** A linear combination of the second derivatives at the left
c ****** and right ends of the last cell is specified.
c
        a(n)=c1(1)
        c(n)=c1(2)
        fpp(n)=c1(3)
c
      case (8)
c
c ****** A linear combination of the first derivatives at the left
c ****** and right ends of the last cell is specified.
c
        dx=x(n)-x(n-1)
        a(n)=dx*(c1(1)/3-c1(2)/6)
        c(n)=dx*(c1(1)/6-c1(2)/3)
        fpp(n)=c1(3)-(c1(1)+c1(2))*(f(n)-f(n-1))/dx
c
      case default
c
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(N).'
        write (*,*) '### IBC1 is invalid.'
        write (*,*) 'IBC1 = ',ibc1
        call exit (1)
c
      end select
c
c-----------------------------------------------------------------------
c ****** Solve the tridiagonal system for the second derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      call trid (n,c,a,b,fpp,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine spline_periodic_type1 (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients for a periodic function.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** This routine assumes that the data in F is periodic, such
c ****** that F(N) = F(1), so that the first and last point are
c ****** repeated and represent the same location.  Thus the
c ****** periodicity length is X(N)-X(1).
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n-1) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 3 points.
c
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the periodic boundary condition at X(1).
c-----------------------------------------------------------------------
c
      dxp=x(2)-x(1)
      dxm=x(n)-x(n-1)
      dxh=dxp+dxm
      a(1)=dxh*third
      c(1)=dxm*sixth
      b(1)=dxp*sixth
      fpp(1)=(f(2)-f(1))/dxp-(f(n)-f(n-1))/dxm
c
c-----------------------------------------------------------------------
c ****** Solve the (periodic) tridiagonal system for the second
c ****** derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      m=n-1

      call trid_periodic (m,c,a,b,fpp,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      fpp(n)=fpp(1)
c
      return
      end
c#######################################################################
      subroutine spline_periodic_type2 (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients for a periodic function.
c
c-----------------------------------------------------------------------
c
c ****** Get the coefficients of a cubic spline interpolant to
c ****** the function values F(i) defined at the points X(i),
c ****** i=1,...,N.
c
c ****** The computed coefficients (which are actually the second
c ****** derivatives of F at the mesh points) are returned in the
c ****** array FPP.
c
c ****** Use routine SPLINT to evaluate the spline at a
c ****** particular position.
c
c-----------------------------------------------------------------------
c
c ****** This routine assumes that the data in F is periodic, such
c ****** that F(N-1) = F(1) and F(N)  = F(2), so that two sets of
c ****** points are repeated and represent the same locations.
c ****** Thus the periodicity length is X(N-1)-X(1).
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(2:n-1) :: a,b,c
c
c-----------------------------------------------------------------------
c
c ****** Check that there are at least 4 points.
c
      if (n.lt.4) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 4 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
c
c ****** Check that the mesh is monotonic.
c
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
c
c-----------------------------------------------------------------------
c ****** Set the coefficients for the tridiagonal solve at the
c ****** internal points.
c-----------------------------------------------------------------------
c
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
c
c-----------------------------------------------------------------------
c ****** Solve the (periodic) tridiagonal system for the second
c ****** derivative.
c-----------------------------------------------------------------------
c
c ****** The second derivative is returned in FPP.
c
      m=n-2

      call trid_periodic (m,c,a,b,fpp(2),ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### The tridiagonal matrix relating the'//
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
c
      fpp(n)=fpp(2)
      fpp(1)=fpp(n-1)
c
      return
      end
c#######################################################################
      subroutine ezspline (n,x,f,fpp)
c
c-----------------------------------------------------------------------
c
c ****** Calculate cubic spline coefficients.
c
c ****** Easy-to-use version of SPLINE.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls SPLINE with boundary conditions
c ****** IBC0=1 and IBC1=1.  See the comments in SPLINE to
c ****** see what this means.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
c
      intent(in) :: n,x,f
      intent(out) :: fpp
c
c-----------------------------------------------------------------------
c
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
c
c-----------------------------------------------------------------------
c
      ibc0=1
      ibc1=1
c
      c0(:)=0.
      c1(:)=0.
c
      call spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
c
      return
      end
c#######################################################################
      function splint (n,x,f,fpp,xv,tab)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a 1D cubic spline.
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the function values F(i) and the second
c ****** derivatives FPP(i), defined at the points X(i),
c ****** i=1,...,N, need to be specified.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval that contains XV.
c
c ****** The routine SPLINE can be used to compute FPP.
c
c ****** This routine uses routine SPEVAL to evaluate the spline.
c
c ****** The value of the spline at XV is returned as the
c ****** function value.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
      use locate_interval_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
      real(r_typ) :: xv
      type(itab), optional :: tab
      real(r_typ) :: splint
c
      intent(in) :: n,x,f,fpp,xv,tab
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: speval
c
c-----------------------------------------------------------------------
c
c ****** Find the mesh interval that encloses XV.
c
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
c
c ****** Evaluate the cubic spline.
c
      splint=speval(x(i),x(i+1),f(i),f(i+1),fpp(i),fpp(i+1),xv)
c
      return
      end
c#######################################################################
      function speval (x1,x2,f1,f2,fpp1,fpp2,xv)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the function values F1 and F2 and the second
c ****** derivatives FPP1 and FPP2, defined at the left and right
c ****** ends of the interval X1 and X2 need to be specified.
c
c ****** The value of the spline at XV is returned as the
c ****** function value.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x1,x2,f1,f2,fpp1,fpp2
      real(r_typ) :: xv
      real(r_typ) :: speval
c
      intent(in) :: x1,x2,f1,f2,fpp1,fpp2,xv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dx,a,b,c,d
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the cubic spline.
c
      dx=x2-x1
c
      b=(xv-x1)/dx
      a=one-b
c
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
c
      speval=a*f1+b*f2+c*fpp1+d*fpp2
c
      return
      end
c#######################################################################
      subroutine speval_get_abcd (x1,x2,xv,a,b,c,d)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c ****** This version splits the calculation into two parts for
c ****** efficiency.
c
c ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
c
c ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
c ****** when multiple spline evaluations are being performed for
c ****** the same position.
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline interpolant at X=XV.
c
c ****** On input, the left and right limits of the spline interval
c ****** X1 and X2, as well as the position at which the spline is
c ****** being evaluated, XV, need to be specified.
c
c ****** The coefficients of the spline A, B, C, and D are
c ****** returned on output.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x1,x2,xv,a,b,c,d
c
      intent(in) :: x1,x2,xv
      intent(out) :: a,b,c,d
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dx
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the coefficients of the cubic spline.
c
      dx=x2-x1
c
      b=(xv-x1)/dx
      a=one-b
c
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
c
      return
      end
c#######################################################################
      function speval_abcd (a,b,c,d,f1,f2,fpp1,fpp2)
c
c-----------------------------------------------------------------------
c
c ****** Evaluate a cubic spline.
c
c ****** This version splits the calculation into two parts for
c ****** efficiency.
c
c ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
c
c ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
c ****** when multiple spline evaluations are being performed for
c ****** the same position.
c
c-----------------------------------------------------------------------
c
c ****** On input, the coefficients A, B, C, and D, and the
c ****** function values F1 and F2 and the second
c ****** derivatives FPP1 and FPP2, need to be specified.
c
c ****** The value of the spline is returned as the
c ****** function value.
c
c ****** The coefficients A, B, C, and D can be obtained using
c ****** routine SPEVAL_GET_ABCD.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: a,b,c,d,f1,f2,fpp1,fpp2
      real(r_typ) :: speval_abcd
c
      intent(in) :: a,b,c,d,f1,f2,fpp1,fpp2
c
c-----------------------------------------------------------------------
c
c ****** Evaluate the cubic spline.
c
      speval_abcd=a*f1+b*f2+c*fpp1+d*fpp2
c
      return
      end
c#######################################################################
      function locate_interval (n,x,xv,tab,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Locate a mesh interval.
c
c-----------------------------------------------------------------------
c
c ****** Find the interval I in table X(i), i=1,2,...,N,
c ****** that encloses the value XV, i.e., such that
c ****** X(I).le.XV.le.X(I+1).
c
c ****** For the special case when N=1, XV must equal X(1)
c ****** exactly, otherwise an error occurs.
c
c ****** The optional argument TAB is an inverse interpolation
c ****** table that can be used to speed up the search for the
c ****** interval.
c
c ****** If the optional argument IERR is specified, then this
c ****** routine will return when an error occurs with IERR=1.
c ****** If no error occurs, IERR=0 is returned.  When IERR is not
c ****** specified, this routine will terminate the program
c ****** with a printed error message.
c
c ****** The mesh interval I is returned as the function value.
c
c-----------------------------------------------------------------------
c
      use number_types
      use invint_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      type(itab), optional :: tab
      integer, optional :: ierr
      integer :: locate_interval
c
      intent(in) :: n,x,xv,tab
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,ig
      real(r_typ) :: xi,fiv,alpha
c
c-----------------------------------------------------------------------
c
      if (present(ierr)) then
        ierr=0
      end if
c
c ****** For the special case when the table has only one
c ****** point (N=1), the inverse table is not used.  In this
c ****** case it is necessary for XV to equal X(1) exactly,
c ****** otherwise this routine exits with an error.
c
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          locate_interval=i
          return
        else
          go to 900
        end if
      end if
c
c ****** Search for the interval depending on whether the optional
c ****** inverse interpolation table TAB was specified.
c
      if (.not.present(tab)) then
c
c ****** Search without an inverse interpolation table.
c
        do i=1,n-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            locate_interval=i
            return
          end if
        enddo
c
      else
c
c ****** Search with an inverse interpolation table.
c
c ****** Get an estimate of the nearest grid point location in
c ****** the (uniform) inverse interpolation table.
c
        xi=one+(xv-x(1))*tab%d
        i=xi
        i=max(i,1)
        i=min(i,tab%n-1)
        alpha=xi-i
        fiv=(one-alpha)*tab%f(i)+alpha*tab%f(i+1)
c
c ****** Set IG to be the guess for the nearest grid point.
c
        ig=fiv
        ig=max(ig,1)
        ig=min(ig,n-1)
c
        if (xv.ge.x(ig)) then
c
c ****** Search forwards.
c
          do i=ig,n-1
            if (xv.ge.x(i).and.xv.le.x(i+1)) then
              locate_interval=i
              return
            end if
          enddo
c
        else
c
c ****** Search backwards.
c
          do i=ig-1,1,-1
            if (xv.ge.x(i).and.xv.le.x(i+1)) then
              locate_interval=i
              return
            end if
          enddo
c
        end if
c
      end if
c
  900 continue
c
c ****** Value not found.
c
c ****** If IERR was passed, set IERR=1 and return; otherwise,
c ****** write an error message and terminate the program.
c
      if (present(ierr)) then
        ierr=1
        return
      else
        write (*,*)
        write (*,*) '### ERROR in LOCATE_INTERVAL:'
        write (*,*) '### The value requested was not found in'//
     &              ' the table:'
        write (*,*) 'Value requested = ',xv
        write (*,*) 'Minimum table value = ',x(1)
        write (*,*) 'Maximum table value = ',x(n)
        call exit (1)
      end if
c
      end
c#######################################################################
      subroutine trid (n,c,a,b,d,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the tridiagonal system of equations:
c
c         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
c
c        for i=2,...,N-1, with
c
c           A(1)*X(1) + B(1)*X(2) = D(1)
c
c        and
c
c           C(N)*X(N-1) + A(N)*X(N) = D(N)
c
c ****** Note that C(1) and B(N) are not referenced.
c
c ****** D is overwritten with the solution.
c
c ****** Return IERR=0 for a successful completion.  If the
c ****** matrix is singular, this routine returns IERR=1 and
c ****** D is invalid.
c
c ****** This routine does not do any pivoting, so the solution
c ****** is not guaranteed to be accurate unless the matrix is
c ****** well-conditioned (e.g., diagonally dominant).
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
c
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
c
c-----------------------------------------------------------------------
c
      ierr=1
c
c ****** Copy A to AA, since it will be overwritten during the
c ****** elimination.  This prevents A from being overwritten.
c
      aa=a
c
c ****** Forward elimination.
c
      if (aa(1).eq.0.) return
      d(1)=d(1)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        d(i)=ace*(d(i)-c(i)*d(i-1))
      enddo
c
c ****** Backward substitution.
c
      do i=n-1,1,-1
        d(i)=d(i)-aa(i)*d(i+1)
      enddo
c
c ****** Set the error return flag to indicate successful completion.
c
      ierr=0
c
      return
      end
c#######################################################################
      subroutine trid_periodic (n,c,a,b,d,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the tridiagonal system of equations:
c
c         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
c
c        for i=2,...,N-1, with
c
c           C(1)*X(N) + A(1)*X(1) + B(1)*X(2) = D(1)
c
c        and
c
c           C(N)*X(N-1) + A(N)*X(N) + B(N)*X(1) = D(N)
c
c ****** D is overwritten with the solution.
c
c ****** Return IERR=0 for a successful completion.  If the
c ****** matrix is singular, this routine returns IERR=1 and
c ****** D is invalid.
c
c ****** This routine does not do any pivoting, so the solution
c ****** is not guaranteed to be accurate unless the matrix is
c ****** well-conditioned (e.g., diagonally dominant).
c
c ****** This system arises for periodic solutions.
c
c ****** This routine uses the Sherman-Morrison formula for
c ****** updating a matrix inverse with a low-rank modification.
c ****** The modification arises from the changes introduced by the
c ****** periodic boundary conditions.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
c
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
      real(r_typ), dimension(n) :: y
      real(r_typ), dimension(n,2) :: z2
      real(r_typ), dimension(2,2) :: t,tinv
      real(r_typ) :: detinv
      real(r_typ), dimension(2) :: tvty
c
c-----------------------------------------------------------------------
c
      ierr=1
c
c ****** First, solve the (non-periodic) system with an RHS
c ****** equal to D.
c
      y=d
c
      call trid (n,c,a,b,y,ierr)
c
      if (ierr.ne.0) return
c
c ****** Next, solve the two systems for the "inhomogenous part".
c
      z2=0.
      z2(1,1)=one
      z2(n,2)=one
c
c ****** Copy A to AA, since it will be overwritten during the
c ****** elimination.  This prevents A from being overwritten.
c
      aa=a
c
c ****** Forward elimination.
c
      if (aa(1).eq.0.) return
      z2(1,:)=z2(1,:)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        z2(i,:)=ace*(z2(i,:)-c(i)*z2(i-1,:))
      enddo
c
c ****** Backward substitution.
c
      do i=n-1,1,-1
        z2(i,:)=z2(i,:)-aa(i)*z2(i+1,:)
      enddo
c
c ****** Invert the 2 x 2 system.
c
      t(1,1)=one+z2(n,1)*b(n)
      t(1,2)=z2(n,2)*b(n)
      t(2,1)=z2(1,1)*c(1)
      t(2,2)=one+z2(1,2)*c(1)
c
      denom=t(1,1)*t(2,2)-t(1,2)*t(2,1)
      if (denom.eq.0.) return
      detinv=one/denom
c
      tinv(1,1)= detinv*t(2,2)
      tinv(2,2)= detinv*t(1,1)
      tinv(1,2)=-detinv*t(1,2)
      tinv(2,1)=-detinv*t(2,1)
c
c ****** Construct the final periodic solution.
c
      tvty(1)=tinv(1,1)*b(n)*y(n)+tinv(1,2)*c(1)*y(1)
      tvty(2)=tinv(2,1)*b(n)*y(n)+tinv(2,2)*c(1)*y(1)
c
      d(:)=y(:)-z2(:,1)*tvty(1)-z2(:,2)*tvty(2)
c
c ****** Set the error return flag to indicate successful completion.
c
      ierr=0
c
      return
      end
c#######################################################################
      subroutine get_dips_map_3d
c
c-----------------------------------------------------------------------
c
c ****** Compute a 3D dips map.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage for the coronal hole map.
c
      real(r_typ), dimension(nrss,ntss,npss) :: dips
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done, private_dips
c
c-----------------------------------------------------------------------
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) '### Computing a 3D dips map:'
      end if
c
c ****** Check that the coronal hole map output file name is not
c ****** blank, since this does not make sense.
c
      if (dips_map_3d_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_DIPS_MAP_3D:'
        write (*,*) '### A coronal hole map was requested, yet'//
     &              ' the output file name is blank.'
        call exit (1)
      end if


      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
c
      n_total=nrss*ntss*npss
      n_completed=0

c
c$omp parallel do
c$omp& private(i,j,k,xfl0,private_dips)
c$omp& private(nc,diag_step,pct_done)
c$omp& collapse(3)
c$omp& schedule(dynamic,iterations_per_thread)
      do i=1,nrss
        do k=1,npss
          do j=1,ntss
c
c ****** Update the iteration counter for diagnostic
c ****** purposes.
c
            if (verbose.gt.0) then
c$omp critical
              n_completed=n_completed+1
              nc=n_completed
c$omp end critical
            end if
c
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
            private_dips=0.
            call get_dip(xfl0,private_dips)
            dips(i,j,k)=private_dips
c
c ****** Write progress diagnostics if requested.
c
            if (verbose.gt.0) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
c
          enddo
        enddo
      enddo
c$omp end parallel do

c
c ****** Write the coronal hole map.
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Writing the 3D coronal hole map to file: ',
     &              trim(dips_map_3d_output_file)
      end if
c
      call wrhdf_3d (dips_map_3d_output_file,.true.,
     &               nrss,ntss,npss,dips,rss,tss,pss,
     &               hdf32,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_DIPS_MAP_3D:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
c
      return
      end
c#######################################################################
      subroutine get_dip(xfl0,dips)
c
c-----------------------------------------------------------------------
c
c ****** Compute a 3D dips
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
c
      real(r_typ) :: dips
c
c-----------------------------------------------------------------------
c
      integer :: is
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      logical :: dipb,dipf
c
c ****** Field line trace storage buffers.
c
      type(traj) :: xtf, xtb
c
c-----------------------------------------------------------------------
c
c ****** Set the tracing direction to be either along the direction
c ****** of the magnetic field or along the directon of increasing
c ****** radius.
c
      ds%direction_is_along_b=trace_slice_direction_is_along_b
c
      ds_f=ds
      ds_f%direction=1
c
      ds_b=ds
      ds_b%direction=-1
c
c ****** Trace a field line in the forward direction along B.
c
      call allocate_trajectory_buffer (xtf)
      call allocate_trajectory_buffer (xtb)
      call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb,xtf)
c
c ****** Check that the field line reached R0 or R1.
c
      if (ttb) then
        f_trace_reached_boundary=.true.
        f_trace_on_r0=xfl1(1).eq.b%lim0(1)
        f_trace_on_r1=xfl1(1).eq.b%lim1(1)
        f_br_positive=bs1(1).ge.0.
      else
        f_trace_reached_boundary=.false.
      end if
c
c ****** Trace a field line in the backward direction along B.
c
      call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb,xtb)
c
c ****** Check that the field line reached R0 or R1.
c
      if (ttb) then
        b_trace_reached_boundary=.true.
        b_trace_on_r0=xfl1(1).eq.b%lim0(1)
        b_trace_on_r1=xfl1(1).eq.b%lim1(1)
        b_br_positive=bs1(1).ge.0.
      else
        b_trace_reached_boundary=.false.
      end if
c
      if (f_trace_reached_boundary.and.
     &      b_trace_reached_boundary) then
        if (f_trace_on_r0.and.b_trace_on_r0) then
          dipf=.false.
          dipb=.false.
          do is=1,min(xtf%npts,ns_dips)
            if(xtf%x(1)%f(is).gt.xfl0(1)) then
              dipf=.true.
              exit
            endif
          enddo
          do is=1,min(xtb%npts,ns_dips)
            if(xtb%x(1)%f(is).gt.xfl0(1)) then
              dipb=.true.
              exit
            endif
          enddo
          if (dipb.and.dipf) then
            dips=one
          endif
        end if
      end if
      call deallocate_trajectory_buffer (xtf)
      call deallocate_trajectory_buffer (xtb)
c
      return
      end
c#######################################################################
      subroutine set_up_integration
c
c-----------------------------------------------------------------------
c
c ****** Set up integration along field lines.
c ****** Read scalar field.
c
c-----------------------------------------------------------------------
c
      use number_types
      use vars
      use params
      use integrate_fl
      use files
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      do_integral_along_fl=.true.
c
c ****** Read the scalar field
c
      if (verbose.gt.0) then
        write (*,*)
        write (*,*) 'Reading data file: ',trim(scalar_input_file)
      end if
c
      call rdhdf (scalar_input_file,scalar_field,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SET_UP_INTEGRATION:'
        write (*,*) '### Could not read scalar field.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(scalar_input_file)
        call exit (1)
      end if
c
      if (scalar_field%ndim.ne.3.or..not.scalar_field%scale)  then
        write (*,*)
        write (*,*) '### ERROR in SET_UP_INTEGRATION:'
        write (*,*) '### Invalid or missing scales in scalar file.'
        write (*,*) 'File name: ',trim(scalar_input_file)
        call exit (1)
      end if
c
      call build_inverse_tables (scalar_field,inv_sf)
c
      return
      end
c#######################################################################
      subroutine magnetic_field_function (time,rtp,s,v)
c
c-----------------------------------------------------------------------
c
c ****** Analytic magnetic field function.
c
c-----------------------------------------------------------------------
c
      use number_types
      use magfld_func_def
      use magfld_func_index
      use magfld_func_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: time
      logical :: rtp
      real(r_typ), dimension(3) :: s
      real(r_typ), dimension(3) :: v
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x,y,z
      real(r_typ) :: r,t,p
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: br_pfss_bkg
      real(r_typ), external :: bt_pfss_bkg
      real(r_typ), external :: bp_pfss_bkg
c
c-----------------------------------------------------------------------
c
      if (rtp) then
        r=s(1)
        t=s(2)
        p=s(3)
      else
        x=s(1)
        y=s(2)
        z=s(3)
      end if
c
      select case (function_index)
      case (FUNC_TYPE_DIPOLE)
        v(1)=two*b0*cos(t)/r**3
        v(2)=b0*sin(t)/r**3
        v(3)=0.
      case (FUNC_TYPE_PFSS_BKG)
        v(1)=br_pfss_bkg(r,t,p,mu,rss)
        v(2)=bt_pfss_bkg(r,t,p,mu,rss)
        v(3)=bp_pfss_bkg(r,t,p,mu,rss)
      case default
        write (*,*)
        write (*,*) '### ERROR in MAGNETIC_FIELD_FUNCTION:'
        write (*,*) '### Invalid function requested:'
        write (*,*) 'FUNCTION_INDEX = ',function_index
        call exit (1)
      end select
c
      return
      end
c#######################################################################
      function br_pfss_bkg (r, theta, phi, mu, Rss)
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r
      real(r_typ) :: theta
      real(r_typ) :: phi
      real(r_typ) :: mu
      real(r_typ) :: Rss
      real(r_typ) :: br_pfss_bkg
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: t1
      real(r_typ) :: t4
      real(r_typ) :: t10
c
c-----------------------------------------------------------------------
c
      t1 = Rss ** 2
      t4 = r ** 2
      t10 = cos(theta)
      br_pfss_bkg = (0.1D1 / t1 / Rss + 0.2D1 / t4 / r) * mu * t10
c
      return
      end
c#######################################################################
      function bt_pfss_bkg (r, theta, phi, mu, Rss)
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r
      real(r_typ) :: theta
      real(r_typ) :: phi
      real(r_typ) :: mu
      real(r_typ) :: Rss
      real(r_typ) :: bt_pfss_bkg
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: t1
      real(r_typ) :: t4
      real(r_typ) :: t9
c
c-----------------------------------------------------------------------
c
      t1 = Rss ** 2
      t4 = r ** 2
      t9 = sin(theta)
      bt_pfss_bkg = (-0.1D1 / t1 / Rss + 0.1D1 / t4 / r) * mu * t9
c
      return
      end
c#######################################################################
      function bp_pfss_bkg (r, theta, phi, mu, Rss)
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r
      real(r_typ) :: theta
      real(r_typ) :: phi
      real(r_typ) :: mu
      real(r_typ) :: Rss
      real(r_typ) :: bp_pfss_bkg
c
c-----------------------------------------------------------------------
c
      bp_pfss_bkg=0.
c
      return
      end
c#######################################################################
c
c Changelog:
c
c        08/30/2005, ZM, Version 1.06:
c
c         - Converted MAPFL v1.05 into a tool.
c         - Improved the tracing accuracy of short field lines
c           (e.g., those near the neutral line).
c         - Improved the field line integrator.
c         - Added the ability to perform a 3D mapping of the
c           field lines.
c
c        06/16/2006, ZM, Version 1.07:
c
c         - Fixed the effect of rondoff errors in generating the
c           meshes for the calculation which caused initial points
c           to lie outside the computation domain.
c
c        06/21/2006, ZM, Version 1.08:
c
c         - Allowed the ability of specifying planes for the
c           field line mapping region and the 3D mapping.
c           This extends the capability of the tool to get
c           Q on 2D slices.
c
c        07/07/2006, ZM, Version 1.09:
c
c         - Added the ability to do new POT3D files,
c           new-old POT3D files, and old POT3D files.
c
c        08/22/2006, ZM, Version 1.10:
c
c         - Corrected a bug in the computation of Q for field
c           lines that wrap periodically at the join between
c           phi=0 and phi=2*pi.
c         - Corrected the computation of Q on open field lines
c           to take into account the difference in the radii
c           of the initial and final field line footpoints.
c         - Added the ability to use cubic spline interpolation.
c           This increases the storage requirements substantially,
c           and makes the computation significantly slower.
c
c        09/17/2006, ZM, Version 1.11:
c
c         - Reverted to the Cartesian field line integrator.
c           The accuracy of the spherical integrator was being
c           cast into doubt, though it was not proven to be
c           bad.
c         - Improved the accuracy with which the final point
c           (i.e., the end point at the r boundaries) was
c           being clipped to the boundary.
c         - Changed the computation of Q to be staggered with
c           respect to the mapping quantities.
c         - Fixed the backward mapping routine to compute Q.
c
c        09/29/2006, ZM, Version 1.12:
c
c         - Improved the field line integrator to use a variable
c           step size.  It is now possible to select either a uniform
c           step size for the field line integration, as before,
c           or to use a variable step size based on the local radius
c           of curvature of the field line, and the local mesh size
c           of the B files along the field line trace.
c         - This has changed the format of the input file.
c
c        01/28/2007, ZM, Version 1.13:
c
c         - Added the ability to compute the mapping on 2D slices
c           through the 3D volume.
c         - This has changed the format of the input file.
c         - Cleaned up some formatting of the output.
c
c        02/19/2008, ZM, Version 1.14:
c
c         - Changed the default behavior that terminated the code
c           when more than 100 bad field line traces were found.
c           By default this is now disabled, but can be put back
c           by setting the variable MAX_BAD_FIELDLINES.
c
c        04/06/2009, ZM, Version 1.15:
c
c         - Added the ability to use an analytic function to define
c           the magnetic field.
c         - Corrected bugs involving how the expansion factor and Q
c           were computed for the backward trace.
c         - Performed a cosmetic clean-up of the code.
c
c        04/21/2009, ZM, Version 1.16:
c
c         - Added the ability to use Slava's new formulation
c           to compute Q directly on slices within the domain
c           by tracing a bundle of field lines from points in
c           the domain forward and backward to the boundaries.
c         - Cleaned up the way mapping along "slices" in the 3D
c           volume is implemented.  These "slices" can now be
c           lines (i.e., 1D files), 2D slices, or 3D volumes.
c           These are all defined by reading in rectilinear
c           files (1D, 2D, or 3D) that contain the (r,t,p) or
c           (x,y,z) starting coordinates for the mapping.
c         - Added a check to make sure that the input file
c           has the correct number of lines.
c
c        02/18/2010, ZM, Version 1.17:
c
c         - Made the ability to trace from a slice more flexible.
c           It is now possible to map from a slice along the
c           forward and backward directions separately.  It is
c           also possible to specify this direction to be either
c           along B or along the direction of increasing radius.
c         - Added the ability to directly compute coronal hole
c           maps at a given radius.  These are coded by magnetic
c           field polarity.
c
c        04/12/2010, ZM, Version 1.18:
c
c         - Added the ability of specifying the r, t, and p meshes
c           to be used for the calculation.  These can be specified
c           as 1D HDF files or as uniform meshes.
c         - Allowed the phi coordinate to be outside the [0,2*pi]
c           range in the input file.  It is now properly wrapped
c           into the main interval during the calculation.
c
c        04/26/2010, ZM, Version 1.19:
c
c         - Added the ability to compute 3D coronal hole maps.
c           These are useful to compute 3D open/closed field regions
c           for use, perhaps, in developing heating masks.
c         - Removed the writing of warning messages about field
c           lines that did not reach the boundaries in the coronal
c           hole map traces.  This was not really necessary since
c           such field lines are already flagged (with values of
c           "-2") in the coronal hole maps.
c
c        02/07/2011, ZM, Version 1.20:
c
c         - Added a multi-threading capability using OpenMP.
c           This version runs in parallel.  This required a
c           restructuring of the code to improve the modularity,
c           to make the code thread-safe, and to improve the
c           parallel performance.
c         - It is necessary to use "thread-safe" programming
c           in the parallel sections from now on.
c
c        05/18/2011, ZM, Version 1.21:
c
c         - Corrected a bug in the periodic wrapping of the phi
c           coordinate when it was outside the range [0,2*pi].
c           It was not being wrapped to the main periodic interval
c           in the GETB routine.
c
c        08/17/2012, ZM, Version 1.22:
c
c         - Fixed a minor bug that was discovered by compiling
c           with GFORTRAN.
c
c        04/29/2013, ZM, Version 1.23:
c
c         - Changed the interpretation of the value of DEBUG_LEVEL
c           for debugging output.  This was done to optimize
c           Slava's output of field line traces.  The functionality
c           of the program was not otherwise modified.
c
c        05/13/2013, ZM, Version 1.24:
c
c         - Added the ability to write field line traces to HDF
c           files.  This option can be selected when tracing from
c           a "slice".  It is not intended to be used when doing
c           very high resolution mapping computations, since it
c           can produce a large amount of data.  To get the HDF
c           files, set DEBUG_LEVEL.ge.2.  The field lines traces
c           for the forward and/or backward trace will be written
c           to HDF files named:
c
c             field_line_traces_forward_rtp.hdf
c             field_line_traces_forward_xyz.hdf
c             field_line_traces_backward_rtp.hdf
c             field_line_traces_backward_xyz.hdf
c
c        09/25/2015, ZM, Version 1.25:
c
c         - Changed the way short field lines are treated.
c           Because of strange behavior noticed by Slava in a
c           certain case, I changed the way the step size for
c           "short" field lines was controlled.  It turned out
c           that a "short" field line could become a very long
c           field line when retraced!  Previously, these field
c           lines were retraced with the miniumum step size, which
c           led to a very long field line with lots of points.
c           I relaxed this requirement on the step size once the
c           number of points in the retraced field line exceeded
c           the minimum number of points.
c
c        01/20/2016, RC, Version 1.26:
c
c         - Modified tracefl() routine to eliminate the use of "goto"s.
c         - Minor speed improvements (~3%).
c         - Manually added spline library into code to allow
c           for the development of future optimizations.
c
c        02/09/2016, ZM, Version 1.26ZM:
c
c         - Fixed the check that adds a phi point to the B arrays
c           for cases when all the B components are on the same
c           mesh.  The code now checks if the point at phi=2*pi is
c           already present, and only adds a point if it is not.
c         - Fixed the storage of field line traces in the field line
c           buffers.  Previously, when "short" field lines were
c           detected, the buffers were not reinitialized when these
c           short field lines were retraced with a smaller step size,
c           leading to incorrect saved field line traces.
c         - Fixed the computation of Q on slices for field line
c           footpoints that approach the poles in routine GETQ.
c           The previous differencing was not accurate for field
c           line footpoints near the poles.  The new scheme switches
c           to Cartesian basis vectors in a small neighborhood of
c           the poles.
c
c        03/25/2016, ZM, Version 1.27ZM:
c
c         - Added code to snap the theta and phi limits of the B
c           meshes to the correct values.  Previously, these could
c           be slightly off due to the precision limitations of
c           32-bit HDF files, leading to possible errors in the
c           tracing of a (very small) subset of field lines.
c
c        04/19/2016, ZM, Version 1.28ZM:
c
c         - Added the ability to write the field line length when
c           mapping forward/backward, and also when computing
c           Q on a slice.
c         - Added the ability to write traces of field lines traced
c           from a slice to individual HDF output files.  These are
c           named:
c
c              <root>_f_####.hdf
c              <root>_b_####.hdf
c
c           where <root> is a string specified in the input file,
c           and #### is a four-digit sequence number (0001, 0002,
c           etc.).  The "f" and "b" labels designate forward/backward
c           traces, as requested in the input file.  There is a flag
c           in the input file to select whether Cartesian (x,y,z)
c           coordinates or spherical (r,t,p) coordinates are
c           written.
c         - Removed the ability to write field lines to a single HDF
c           file when the debug level was set to 2.  This is
c           superseded by the new capability to write the traces to
c           individual HDF files.
c         - COMPTIBILITY CHANGE:  The changes in this version have
c           modified the format
c           of the input file, so unfortunately this version is
c           not backward compatible.
c
c        01/17/2018, ZM, Version 1.29ZM:
c
c         - Changed the way the domain limits are set from the
c           scales in the magnetic field files.  This is an attempt
c           to make sure that roundoff does not cause points that
c           are near theta=pi and phi=2*pi to end up outside the
c           domain.
c
c        01/18/2018, RC, Version 1.30:
c
c         - Merged in 1.26ZM, 1.27ZM, 1.28ZM, and 1.29ZM changes.
c         - Removed the small speed improvement from 1.26 for
c           consistency reasons.
c
c        03/19/2018, RC, Version 1.31:
c
c         - Small change to allow HDF5 trace output files if
c           input b files are HDF5.
c
c        05/10/2018, RC, Version 1.32:
c
c         - Added OpenMP to 3D cubic spline coef calculation.
c
c        05/20/2019, RL, Version 1.33:
c
c         - Added 3D dips map.
c
c        06/28/2019, RL, Version 2.00:
c
c         - COMPATIBILITY CHANGE!!!!
c           Replaced formatted input file with namelist!
c
c        07/01/2019, RC, Version 2.01:
c
c         - Added option to output Slava's signed-log10
c           of Q for forward and backward traces. To use, set
c           SLOGQFFILE and/or SLOGQBFILE to desired output filenames.
c         - Set default values for all namelist parameters and updated
c           mapfl.dat to reflect the default values.
c         - Added error-checking for namelist.
c
c        07/11/2019, RC, Version 2.02:
c
c         - Fixed bug in writing out backward mapping slogq.
c
c        07/12/2019, RC, Version 2.03:
c
c         - Fixed bug in writing out backward mapping slogq (again).
c
c        08/19/2019, RL, Version 2.04:
c
c         - Introduced capability to integrate scalar field along
c           field lines.
c           Speficy INTEGRATE_ALONG_FL=.true. and the name of the
c           file with the field in SCALAR_INPUT_FILE
c           Either TRACE_FWD, or TRACE_BWD, or TRACE_SLICE and
c           COMPUTE_Q_ON_SLICE, must be set true.
c           Specify the file with integrated quantity in LFFILE
c           (for TRACE_FWD or TRACE_BWD true) or in
c           SLICE_LENGTH_OUTPUT_FILE (for TRACE_SLICE and
c           COMPUTE_Q_ON_SLICE true).
c
c        03/02/2021, AP,CD Version 2.0.5:
c
c         - Small modifications for Python/f2py cross compilation.
c         - Debug statements in mesh detection, tweak to OLD MAS check.
c         - Changed version numbering to standard style.
c
c        05/10/2022, EM,CD Version 2.0.6:
c
c         - Changed logic for slice mapping. If compute_q_on_slice
c           is true, then the direct mapping will also be done after.
c         - This makes it easier to get everything at once.
c
c        04/16/2024, RC Version 2.1.0:
c
c         - Moved changelog to bottom fo code.
c         - Changed verbose to a namelist parameter and an integer
c           (set greater than 0 to activate).
c         - Changed stats variable to an integer to match verbose.
c         - Changed command line.  Now, use as: mapfl INFILE
c           If INFILE is not supplied, it defaults to "mapfl.in".
c         - Added number_types module.
c         - Integrated analytic magnetic field modules and routines
c           into main code and namelist.
c           Note that "rss" is no longer a parameter for the feature
c           (it was not being used anyways).
c         - Updated intents and dummy variables to avoid
c           argument type mismatches.
c
c-----------------------------------------------------------------------
c