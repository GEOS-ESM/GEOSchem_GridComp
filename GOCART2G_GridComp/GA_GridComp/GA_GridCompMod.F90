#include "MAPL_Generic.h"

module GA_GridCompMod

   use ESMF
   use MAPL
   use Chem_MieTableMod2G

   implicit none
   private

   public :: GA_GridComp

   type :: GA_GridComp
       type(Chem_Mie), dimension(2)    :: rad_MieTable, diag_MieTable
       real, allocatable      :: radius(:)      ! particle effective radius [um]
       real, allocatable      :: rlow(:)        ! particle effective radius lower bound [um]
       real, allocatable      :: rup(:)         ! particle effective radius upper bound [um]
       real, allocatable      :: rhop(:)        ! soil class density [kg m-3]
       real, allocatable      :: fscav(:)       ! scavenging efficiency
!      logical                :: scav_byColdCloud ! new flag example
       real, allocatable      :: molwght(:)     ! molecular weight            !NOT UNIVERSAL ONLY FOR GASES, 
       real, allocatable      :: fnum(:)        ! number of particles per kg mass
       integer                :: rhFlag
       integer                :: nbins
       integer                :: km             ! vertical grid dimension
       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: instance       ! data or computational instance
    contains
       procedure :: load_from_config
    end type GA_GridComp


    !LOCALS
     integer :: status
     integer :: nbins

 contains


    subroutine load_from_config(self, cfg, rc)
       class(GA_GridComp), intent(inout) :: self
       type(ESMF_Config), intent(inout) :: cfg
       integer, optional, intent(out) :: rc

       !   Get nbins from cfg
       call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
       nbins = self%nbins
       
       !   Parse config file into private internal state
       !   ----------------------------------------------
       allocate(self%radius(nbins), self%rlow(nbins), self%rup(nbins), self%rhop(nbins), &
                self%fscav(nbins), self%molwght(nbins), self%fnum(nbins), &
                __STAT__)
       
       call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rup,        label='radius_upper:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rhop,       label='soil_density:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%fnum,       label='fnum:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%rhFlag,     label='rhFlag:', __RC__)
       
    end subroutine load_from_config
    
   
end module GA_GridCompMod
