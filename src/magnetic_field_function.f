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
