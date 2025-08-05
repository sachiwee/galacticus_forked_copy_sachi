!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.
!+    Contributions to this file made by: Sachi Weerasooriya
  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use    :: Accretion_Disks                     , only : accretionDisksClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorEddingtonRatio">
   <description>
     A node property extractor which extracts a list of all super-massive black hole eddington ratios
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorEddingtonRatio
     !!{
     A property extractor which extracts a list of all super-massive black hole eddington ratios
     !!}
     private
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class(accretionDisksClass), pointer :: accretionDisks_ => null()
   contains
     final     ::                 eddingtonRatioDestructor
     procedure :: elementCount => eddingtonRatioElementCount
     procedure :: extract      => eddingtonRatioExtract
     procedure :: names        => eddingtonRatioNames
     procedure :: descriptions => eddingtonRatioDescriptions
     procedure :: unitsInSI    => eddingtonRatioUnitsInSI
  end type nodePropertyExtractorEddingtonRatio

  interface nodePropertyExtractorEddingtonRatio
     !!{
     Constructors for the \refClass{nodePropertyExtractorEddingtonRatio} output extractor class.
     !!}
     module procedure eddingtonRatioConstructorParameters
     module procedure eddingtonRatioConstructorInternal
  end interface nodePropertyExtractorEddingtonRatio

contains

  function eddingtonRatioConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorEddingtonRatio} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorEddingtonRatio)                :: self
    type (inputParameters                                 ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass                     ), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass                             ), pointer        :: accretionDisks_          
    !![
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks" name="accretionDisks_" source="parameters"/>
    !!]
    self=nodePropertyExtractorEddingtonRatio(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"/>
    !!]
    return
  end function eddingtonRatioConstructorParameters

  function eddingtonRatioConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorEddingtonRatio} node operator class.
    !!}
    implicit none
    type (nodePropertyExtractorEddingtonRatio)                        :: self
    class(blackHoleAccretionRateClass                     ), intent(in   ), target :: blackHoleAccretionRate_
    class(accretionDisksClass                             ), intent(in   ), target :: accretionDisks_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_,*accretionDisks_"/>
    !!]

    return
  end function eddingtonRatioConstructorInternal

  subroutine eddingtonRatioDestructor(self)
    !!{
    Destructor for the critical overdensity eddingtonRatio set barrier class.
    !!}
    implicit none
    type(nodePropertyExtractorEddingtonRatio), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine eddingtonRatioDestructor

  integer function eddingtonRatioElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorEddingtonRatio), intent(inout) :: self

    eddingtonRatioElementCount=1
    return
  end function eddingtonRatioElementCount

  function eddingtonRatioExtract(self,node,instance) result(eddingtonRatio)
    !!{
    Implement an output extractor for the radiative efficiencies of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    use            :: Numerical_Constants_Math            , only : Pi
    use            :: Numerical_Constants_Atomic          , only : atomicMassUnit
    use            :: Numerical_Constants_Physical        , only : gravitationalConstant, speedLight, thomsonCrossSection
    use            :: Numerical_Constants_Astronomical    , only : massSolar, gigaYear
    implicit none
    double precision                                                  , dimension(:,:), allocatable :: eddingtonRatio
    class           (nodePropertyExtractorEddingtonRatio), intent(inout)               :: self
    type            (treeNode                                        ), intent(inout)               :: node
    type            (multiCounter                                    ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                          )                , pointer     :: blackHole
    integer                                                                                         :: i                                  , countBlackHoles
    double precision                                                                                :: rateMassAccretionSpheroid          , rateMassAccretionHotHalo, &
        &                                                                                              rateMassAccretionNuclearStarCluster,luminosityBolometricAGN ,  &
         &                                                                                             normalization, radiativeEfficiency,                             &
         &                                                                                             rateAccretionBlackHole,bolometricLuminosity, eddingtonLuminosity,  &
         &                                                                                             massProton=1.0072764665789*atomicMassUnit,mass_blackHole
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(eddingtonRatio(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole                =>  node%blackHole(instance=i)
       mass_blackHole      =  blackHole%mass     (          )
       call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNuclearStarCluster)

       ! Get the radiative efficiency of black hole.
       radiativeEfficiency=self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionSpheroid+rateMassAccretionHotHalo+rateMassAccretionNuclearStarCluster)
       
       !Calculate bolometric luminosity
       bolometricLuminosity   =  (speedLight            **2.0d0)          &
         &                          *radiativeEfficiency                  &
         &                          * (  massSolar                        &
         &                         /gigaYear )                            &
         &                         *(                                     &
         &                           +rateMassAccretionHotHalo            &
         &                           +rateMassAccretionSpheroid           &
         &                           +ratemassAccretionNuclearStarCluster &
         &                          )

       ! Calculate Eddington ratio
       eddingtonLuminosity = (4*Pi*gravitationalConstant * massProton* speedLight*mass_blackHole*massSolar) &
         &                     / thomsonCrossSection                                                         
       eddingtonRatio(i,1) =bolometricLuminosity/eddingtonLuminosity
    end do

    return
  end function eddingtonRatioExtract

  subroutine eddingtonRatioNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily eddingtonRatio} properties.
    !!}
    implicit none
    class(nodePropertyExtractorEddingtonRatio), intent(inout)                             :: self
    type (varying_string                                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('eddingtonRatio')
    return
  end subroutine eddingtonRatioNames

  subroutine eddingtonRatioDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily eddingtonRatio} properties.
    !!}
    implicit none
    class(nodePropertyExtractorEddingtonRatio), intent(inout)                             :: self
    type (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Eddington ratios of super-massive black holes in this galaxy.')
    return
  end subroutine eddingtonRatioDescriptions

  function eddingtonRatioUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily eddingtonRatio} properties in the SI system.
    !!}
    implicit none
    double precision                                                  , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorEddingtonRatio), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=1.0d0
    return
  end function eddingtonRatioUnitsInSI
