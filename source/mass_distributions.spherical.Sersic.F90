!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !!{
  Implementation of a S\'ersic mass distribution class.
  !!}

  use    :: Numerical_Interpolation, only : interpolator
  !$ use :: OMP_Lib                , only : omp_lock_kind

  !![
  <massDistribution name="massDistributionSersic">
   <description>A S\'ersic mass distribution class.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionSersic
     !!{
     The S\'ersic density profile.
     !!}
     double precision                                               :: densityNormalization         , mass              , &
          &                                                            radiusHalfMass_              , index_
     ! Tabulation of the Sérsic profile.
     double precision                                               :: coefficient                  , radiusStart
     logical                                                        :: tableInitialized     =.false.
     integer                                                        :: tableCount
     double precision                                               :: tableRadiusMaximum           , tableRadiusMinimum
     double precision                                               :: table3dRadiusHalfMass
     double precision                                               :: table2dRadiusHalfMass
     double precision                   , allocatable, dimension(:) :: tableDensity                 , tableEnclosedMass , &
          &                                                            tablePotential               , tableRadius
     type            (interpolator     )                            :: tableInterpolator
     !$ integer      (omp_lock_kind    )                            :: tableLock
   contains
     !![
     <methods>
       <method description="Tabulate the Sersic profile." method="tabulate" />
       <method description="Return the half mass radius of the profile in projection." method="radiusHalfMassProjected" />
     </methods>
     !!]
     procedure :: tabulate                => sersicTabulate
     procedure :: density                 => sersicDensity
     procedure :: densityRadialMoment     => sersicDensityRadialMoment
     procedure :: massEnclosedBySphere    => sersicMassEnclosedBySphere
     procedure :: massTotal               => sersicMassTotal
     procedure :: potential               => sersicPotential
     procedure :: radiusHalfMass          => sersicRadiusHalfMass
     procedure :: radiusHalfMassProjected => sersicRadiusHalfMassProjected
  end type massDistributionSersic

  interface massDistributionSersic
     !!{
     Constructors for the {\normalfont \ttfamily sersic} mass distribution class.
     !!}
     module procedure sersicConstructorParameters
     module procedure sersicConstructorInternal
  end interface massDistributionSersic

  ! Table granularity for Sersic profiles.
  integer                        , parameter :: tablePointsPerDecade=1000

  ! Module scope variables used in integration and root finding.
  class  (massDistributionSersic), pointer   :: self_
  !$omp threadprivate(self_)

contains

  function sersicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sersic} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSersic)                :: self
    type            (inputParameters       ), intent(inout) :: parameters
    double precision                                        :: mass         , radiusHalfMass, &
         &                                                     index_
    logical                                                 :: dimensionless
    type            (varying_string        )                :: componentType
    type            (varying_string        )                :: massType

    !![
    <inputParameter>
      <name>index</name>
      <variable>index_</variable>
      <defaultValue>4.0d0</defaultValue>
      <description>The S\'ersic index.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusHalfMass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The half mass radius of the S\'ersic profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the S\'ersic profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the S\'ersic profile is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <conditionalCall>
     <call>self=massDistributionSersic(index_,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="radiusHalfMass" value="radiusHalfMass"      parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"          parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sersicConstructorParameters

  function sersicConstructorInternal(index,radiusHalfMass,mass,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``sersic'' mass distribution class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (massDistributionSersic      )                          :: self
    double precision                              , intent(in   )           :: index
    double precision                              , intent(in   ), optional :: radiusHalfMass, mass
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]

    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! Initialize state.
    self%tableInitialized     =.false.
    self%tableRadiusMaximum   =1.0d+3
    self%tableRadiusMinimum   =1.0d-3
    self%table3dRadiusHalfMass=1.0d+0
    self%index_               =index
    !$ call OMP_Init_Lock(self%tableLock)
    ! Tabulate the profile.
    call self%tabulate()
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(radiusHalfMass      )) then
          if (Values_Differ(radiusHalfMass,1.0d0,absTol=1.0d-6)) call Error_Report('radiusHalfMass should be unity for a dimensionless profile (or simply do not specify a half mass radius)'//{introspection:location})
       end if
       if (present(mass                )) then
          if (Values_Differ(mass          ,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'                      //{introspection:location})
       end if
       self%radiusHalfMass_=1.0d0
       self%mass          =1.0d0
    else
       if (present(radiusHalfMass)) then
          self%radiusHalfMass_=radiusHalfMass
       else
          call Error_Report('"radiusHalfMass" must be specified'//{introspection:location})
       end if
       if (present(mass)) then
          self%mass          =mass
       else
          call Error_Report('"mass" must be specified'          //{introspection:location})
       end if
    end if
    return
  end function sersicConstructorInternal

  double precision function sersicDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a S\'ersic mass distribution.
    !!}
    use :: Coordinates , only : assignment(=)               , coordinateSpherical
    implicit none
    class           (massDistributionSersic      ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (coordinateSpherical         )                          :: position
    double precision                                                        :: r

    if (.not.self%matches(componentType,massType)) then
       sersicDensity=0.0d0
       return
    end if
    ! Get position in spherical coordinate system.
    position= coordinates
    ! Compute the density at this position.
    !$ call OMP_Set_Lock(self%tableLock)
    r       =+position%r             () &
         &   /self    %radiusHalfMass_
    call self%tabulate(r)
    sersicDensity=+self%mass                                               &
         &        /self%radiusHalfMass_**3                                 &
         &        *self%tableInterpolator%interpolate(r,self%tableDensity)
    !$ call OMP_Unset_Lock(self%tableLock)
    return
  end function sersicDensity

  double precision function sersicDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Returns a radial density moment for the S\'ersic mass distribution.
    !!}
    implicit none
    class           (massDistributionSersic      ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum          , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    integer                                                                 :: iRadius
    double precision                                                        :: deltaRadius            , integrand              , &
         &                                                                     previousIntegrand      , fractionalRadiusMinimum, &
         &                                                                     fractionalRadiusMaximum

    if (.not.self%matches(componentType,massType)) then
       sersicDensityRadialMoment=0.0d0
       return
    end if
    isInfinite               =.false.
    sersicDensityRadialMoment=0.0d0
    !$ call OMP_Set_Lock(self%tableLock)
    if (present(radiusMinimum)) then
       fractionalRadiusMinimum=radiusMinimum/self%radiusHalfMass_
    else
       fractionalRadiusMinimum=0.0d0
    end if
    if (present(radiusMaximum)) then
       fractionalRadiusMaximum=radiusMaximum/self%radiusHalfMass_
    else
       fractionalRadiusMaximum=self%tableRadius(size(self%tableRadius))
    end if
    do iRadius=1,self%tableCount
       if (iRadius == 1) then
          deltaRadius      =+max(min(self%tableRadius(iRadius  ),fractionalRadiusMaximum),fractionalRadiusMinimum) &
               &            -max(min(                      0.0d0,fractionalRadiusMaximum),fractionalRadiusMinimum)
          previousIntegrand=+0.0d0
       else
          deltaRadius      =+max(min(self%tableRadius(iRadius  ),fractionalRadiusMaximum),fractionalRadiusMinimum) &
               &            -max(min(self%tableRadius(iRadius-1),fractionalRadiusMaximum),fractionalRadiusMinimum)
          previousIntegrand=+self%tableRadius (iRadius-1)**moment &
               &            *self%tableDensity(iRadius-1)
       end if
       integrand                =+self%tableRadius (iRadius)**moment &
            &                    *self%tableDensity(iRadius)
       sersicDensityRadialMoment=+sersicDensityRadialMoment &
            &                    +0.5d0                     &
            &                    *(                         &
            &                      +previousIntegrand       &
            &                      +        integrand       &
            &                     )                         &
            &                    *deltaRadius
    end do
    !$ call OMP_Unset_Lock(self%tableLock)
    sersicDensityRadialMoment=+sersicDensityRadialMoment           &
         &                    *self%mass                           &
         &                    *self%radiusHalfMass_**(moment-3.0d0)
    return
  end function sersicDensityRadialMoment

  double precision function sersicMassTotal(self,componentType,massType)
    !!{
    Computes the total mass for S\'ersic mass distributions.
    !!}
    implicit none
    class(massDistributionSersic      ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       sersicMassTotal=0.0d0
       return
    end if
    sersicMassTotal=self%mass
    return
  end function sersicMassTotal
  
  double precision function sersicMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for S\'ersic mass distributions.
    !!}
    implicit none
    class           (massDistributionSersic      ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: fractionalRadius

    if (.not.self%matches(componentType,massType)) then
       sersicMassEnclosedBySphere=0.0d0
       return
    end if
    if (radius <= 0.0d0) then
       sersicMassEnclosedBySphere=0.0d0
    else
       !$ call OMP_Set_Lock(self%tableLock)
       fractionalRadius=+     radius         &
            &           /self%radiusHalfMass_
       call self%tabulate(fractionalRadius)
       if (fractionalRadius < self%tableRadius(self%tableCount)) then
          sersicMassEnclosedBySphere=+self%mass                                                                   &
               &                     *self%tableInterpolator%interpolate(fractionalRadius,self%tableEnclosedMass)
       else
          sersicMassEnclosedBySphere=+self%mass
       end if
       !$ call OMP_Unset_Lock(self%tableLock)
    end if
    return
  end function sersicMassEnclosedBySphere

  double precision function sersicPotential(self,coordinates,componentType,massType)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a S\'ersic mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)                  , coordinateSpherical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionSersic      ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (coordinateSpherical         )                          :: position
    double precision                                                        :: r

    if (.not.self%matches(componentType,massType)) then
       sersicPotential=0.0d0
       return
    end if
    ! Get position in spherical coordinate system.
    position=coordinates
    ! Compute the potential at this position.
    !$ call OMP_Set_Lock(self%tableLock)
    r       =+position%r             () &
         &   /self    %radiusHalfMass_
    call self%tabulate(r)
    if (r < self%tableRadius(self%tableCount)) then
       sersicPotential=+self%mass                                                 &
            &          /self%radiusHalfMass_                                      &
            &          *self%tableInterpolator%interpolate(r,self%tablePotential)
    else
       sersicPotential=0.0d0
    end if
    !$ call OMP_Unset_Lock(self%tableLock)
    if (.not.self%isDimensionless()) sersicPotential=+gravitationalConstantGalacticus &
         &                                           *sersicPotential
    return
  end function sersicPotential

  double precision function sersicRadiusHalfMass(self,componentType,massType)
    !!{
    Return the half-mass radius of a S\'ersic mass distribution.
    !!}
    implicit none
    class(massDistributionSersic      ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       sersicRadiusHalfMass=0.0d0
       return
    end if
    !$ call OMP_Set_Lock(self%tableLock)
    sersicRadiusHalfMass=+self%radiusHalfMass_
    !$ call OMP_Unset_Lock(self%tableLock)
    return
  end function sersicRadiusHalfMass

  double precision function sersicRadiusHalfMassProjected(self,componentType,massType)
    !!{
    Return the half-mass radius in projection of a S\'ersic mass distribution.
    !!}
    implicit none
    class(massDistributionSersic      ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       sersicRadiusHalfMassProjected=0.0d0
       return
    end if
    sersicRadiusHalfMassProjected=+self%radiusHalfMass_        &
         &                        *self%table2dRadiusHalfMass
    return
  end function sersicRadiusHalfMassProjected

  subroutine sersicTabulate(self,radius)
    !!{
    Tabulate the density and enclosed mass in a dimensionless S\'ersic profile.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range                  , rangeTypeLogarithmic
    use :: Root_Finder             , only : rangeExpandMultiplicative   , rootFinder
    use :: Table_Labels            , only : extrapolationTypeExtrapolate
    implicit none
    class           (massDistributionSersic), intent(inout), target   :: self
    double precision                        , intent(in   ), optional :: radius
    double precision                        , parameter               :: radiusMaximumToTabulate=1.0d+3
    double precision                        , parameter               :: coefficientTolerance   =1.0d-6
    double precision                        , parameter               :: coefficientGuess       =7.67d0
    type            (rootFinder            ), save                    :: finder
    logical                                 , save                    :: finderConstructed      =.false.
    !$omp threadprivate(finder,finderConstructed)
    type            (interpolator          ), allocatable             :: interpolator_
    logical                                                           :: rebuildTable                   , tableHasSufficientExtent
    integer                                                           :: iRadius
    double precision                                                  :: deltaRadius                    , integrand                , &
         &                                                               massPrevious                   , previousIntegrand        , &
         &                                                               radiusActual                   , radiusInfinity
    type            (integrator            )                          :: integrator_

    ! Check if a radius was specified. Use it if so, otherwise use a midpoint radius.
    if (present(radius)) then
       radiusActual=min(radius,radiusMaximumToTabulate)
    else
       radiusActual=sqrt(self%tableRadiusMinimum*self%tableRadiusMaximum)
    end if
    ! Determine if the table must be rebuilt.
    if (self%tableInitialized) then
       rebuildTable=(radiusActual*self%table3dRadiusHalfMass < self%tableRadiusMinimum) &
            &        .or.                                                               &
            &       (radiusActual*self%table3dRadiusHalfMass > self%tableRadiusMaximum)
    else
       rebuildTable=.true.
    end if
    ! Rebuild the table if necessary.
    if (rebuildTable) then
       ! Set a module-scope pointer to self.
       self_ => self
       ! Initialize our root finder.
       if (.not.finderConstructed) then
          finder           =rootFinder(                                               &
               &                       rootFunction       =sersicCoefficientRoot    , &
               &                       toleranceAbsolute  =0.0d0                    , &
               &                       toleranceRelative  =coefficientTolerance     , &
               &                       rangeExpandDownward=0.5d0                    , &
               &                       rangeExpandUpward  =2.0d0                    , &
               &                       rangeExpandType    =rangeExpandMultiplicative  &
               &                      )
          finderConstructed=.true.
       end if
       ! Try building the table until it has sufficient extent to encompass the requested radius.
       tableHasSufficientExtent=.false.
       do while (.not.tableHasSufficientExtent)
          ! Find suitable radius limits.
          self%tableRadiusMinimum=min(self%tableRadiusMinimum,0.5d0*radiusActual*self%table3dRadiusHalfMass)
          self%tableRadiusMaximum=max(self%tableRadiusMaximum,2.0d0*radiusActual*self%table3dRadiusHalfMass)
          ! Determine the number of points at which to tabulate the profile.
          self%tableCount=int(log10(self%tableRadiusMaximum/self%tableRadiusMinimum)*dble(tablePointsPerDecade))+1
          ! Allocate arrays for storing the tables.
          if (allocated(self%tableRadius)) then
             deallocate(self%tableRadius      )
             deallocate(self%tableDensity     )
             deallocate(self%tableEnclosedMass)
             deallocate(self%tablePotential   )
          end if
          allocate(self%tableRadius      (self%tableCount))
          allocate(self%tableDensity     (self%tableCount))
          allocate(self%tableEnclosedMass(self%tableCount))
          allocate(self%tablePotential   (self%tableCount))
          ! Create an array of logarithmically distributed radii.
          self%tableRadius=Make_Range(self%tableRadiusMinimum,self%tableRadiusMaximum,self%tableCount,rangeType=rangeTypeLogarithmic)
          ! Compute the coefficient appearing in the Sérsic profile.
          self%coefficient=finder%find(rootGuess=coefficientGuess)
          ! Compute a suitably large approximation to infinite radius for use in integration.
          radiusInfinity=10.0d0*self%tableRadiusMaximum
          ! Loop over radii and compute the inverse Abel integral required to get the 3D Sérsic profile.
          integrator_=integrator(sersicAbelIntegrand,toleranceRelative=1.0d-3)
          do iRadius=1,self%tableCount
             self%radiusStart          =self       %tableRadius(iRadius)
             self%tableDensity(iRadius)=integrator_%integrate  (                     &
                  &                                             self%radiusStart   , &
                  &                                                  radiusInfinity  &
                  &                                            )
             ! Accumulate the enclosed mass using a simple trapezoidal integration.
             if (iRadius == 1) then
                deltaRadius      =self%tableRadius(iRadius)
                massPrevious     =0.0d0
                previousIntegrand=0.0d0
             else
                deltaRadius      =+self%tableRadius      (iRadius  ) &
                     &            -self%tableRadius      (iRadius-1)
                massPrevious     =+self%tableEnclosedMass(iRadius-1)
                previousIntegrand=4.0d0*Pi*(self%tableRadius(iRadius-1)**2)*self%tableDensity(iRadius-1)
             end if
             integrand           =4.0d0*Pi*(self%tableRadius(iRadius  )**2)*self%tableDensity(iRadius  )
             self%tableEnclosedMass(iRadius)=0.5d0*(previousIntegrand+integrand)*deltaRadius+massPrevious
          end do
          ! Normalize the mass and density to unit total mass.
          self%tableDensity     =self%tableDensity     /self%tableEnclosedMass(self%tableCount)
          self%tableEnclosedMass=self%tableEnclosedMass/self%tableEnclosedMass(self%tableCount)
          ! Find the half mass radius.
          allocate(interpolator_)
          interpolator_             =interpolator(                                                                &
               &                                  self%tableEnclosedMass(1:maxloc(self%tableEnclosedMass,dim=1)), &
               &                                  self%tableRadius      (1:maxloc(self%tableEnclosedMass,dim=1))  &
               &                                 )
          self%table3dRadiusHalfMass=interpolator_%interpolate(0.5d0)
          deallocate(interpolator_)
          ! Scale radii and densities to be in units of the 3D half mass radius.
          self%tableRadius =self%tableRadius /self%table3dRadiusHalfMass
          self%tableDensity=self%tableDensity*self%table3dRadiusHalfMass**3
          ! Store the 2d half mass radius.
          self%table2dRadiusHalfMass=1.0d0/self%table3dRadiusHalfMass
          ! Compute the gravitational potential at each radius using a simple trapezoidal rule integration.
          self%tablePotential(self%tableCount)=0.0d0 ! Assume zero potential at effective infinity.
          do iRadius=self%tableCount-1,1,-1
             self%tablePotential(iRadius)= self%tablePotential(iRadius+1)                                      &
                  &                        -(                                                                  &
                  &                           self%tableEnclosedMass(iRadius+1)/self%tableRadius(iRadius+1)**2 &
                  &                          +self%tableEnclosedMass(iRadius  )/self%tableRadius(iRadius  )**2 &
                  &                         )                                                                  &
                  &                        *0.5d0                                                              &
                  &                        *(                                                                  &
                  &                           self%tableRadius      (iRadius+1)                                &
                  &                          -self%tableRadius      (iRadius  )                                &
                  &                         )
          end do
          ! Test that the table has sufficient extent for the requested radius.
          tableHasSufficientExtent=(radiusActual >= self%tableRadiusMinimum) &
               &                    .and.                                    &
               &                   (radiusActual <= self%tableRadiusMaximum)
       end do
       ! Build the interpolator.
       self%tableInterpolator=interpolator(self%tableRadius,extrapolationType=extrapolationTypeExtrapolate)
       ! Flag that the table is initialized.
       self%tableInitialized=.true.
    end if
    return
  end subroutine sersicTabulate

  double precision function sersicCoefficientRoot(coefficient)
    !!{
    Root function used in finding the coefficient for S\'ersic profiles.
    !!}
    use :: Gamma_Functions, only : Gamma_Function_Incomplete
    implicit none
    double precision, intent(in   ) :: coefficient

    sersicCoefficientRoot=Gamma_Function_Incomplete(2.0d0*self_%index_,coefficient)-0.5d0
    return
  end function sersicCoefficientRoot

  double precision function sersicAbelIntegrand(radius)
    !!{
    The integrand in the Abel integral used to invert the S\'ersic profile to get the corresponding 3-D profile.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    if (radius > self_%radiusStart) then
       sersicAbelIntegrand=+      self_%coefficient*(radius**(1.0d0/dble(self_%index_) -1.0d0)) &
            &              * exp(-self_%coefficient*(radius**(1.0d0/dble(self_%index_))-1.0d0)) &
            &              /                                        dble(self_%index_)          &
            &              /sqrt(                                                               &
            &                    +     radius     **2                                           &
            &                    -self_%radiusStart**2                                          &
            &                   )                                                               &
            &              /Pi
    else
       sersicAbelIntegrand=0.0d0
    end if
    return
  end function sersicAbelIntegrand
