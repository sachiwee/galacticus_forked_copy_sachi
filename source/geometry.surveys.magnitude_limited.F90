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

!!{
Implements a simple full-sky, magnitude limited survey geometry.
!!}

  !![
  <surveyGeometry name="surveyGeometryMagnitudeLimited">
   <description>Implements the geometry of the SDSS survey of \cite{montero-dorta_sdss_2009}.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryFullSky) :: surveyGeometryMagnitudeLimited
     private
     double precision :: magnitudeApparentMinimum, magnitudeApparentMaximum
   contains
     procedure :: distanceMinimum => magnitudeLimitedDistanceMinimum
     procedure :: distanceMaximum => magnitudeLimitedDistanceMaximum
  end type surveyGeometryMagnitudeLimited

  interface surveyGeometryMagnitudeLimited
     !!{
     Constructors for the \refClass{surveyGeometryMagnitudeLimited} survey geometry class.
     !!}
     module procedure magnitudeLimitedConstructorParameters
     module procedure magnitudeLimitedConstructorInternal
  end interface surveyGeometryMagnitudeLimited

contains

  function magnitudeLimitedConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{surveyGeometryMagnitudeLimited} survey geometry class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryMagnitudeLimited)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    double precision                                                :: redshiftMinimum         , redshiftMaximum         , &
         &                                                             magnitudeApparentMinimum, magnitudeApparentMaximum

    !![
    <inputParameter>
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <defaultValue>0.00208d0</defaultValue>
      <description>The minimum redshift for the survey.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>huge(0.38582d0)</defaultValue>
      <source>parameters</source>
      <description>The maximum redshift for the survey.</description>
    </inputParameter>
    <inputParameter>
      <name>magnitudeApparentMinimum</name>
      <defaultValue>14.5d0</defaultValue>
      <source>parameters</source>
      <description>The minimum apparent magnitude for the survey.</description>
    </inputParameter>
    <inputParameter>
      <name>magnitudeApparentMaximum</name>
      <defaultValue>17.77d0</defaultValue>
      <source>parameters</source>
      <description>The maximum apparent magnitude for the survey.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=surveyGeometryMagnitudeLimited(redshiftMinimum,redshiftMaximum,magnitudeApparentMinimum,magnitudeApparentMaximum,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function magnitudeLimitedConstructorParameters

  function magnitudeLimitedConstructorInternal(redshiftMinimum,redshiftMaximum,magnitudeApparentMinimum,magnitudeApparentMaximum,cosmologyFunctions_) result (self)
    !!{
    Default constructor for the \refClass{surveyGeometryMagnitudeLimited} survey geometry class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (surveyGeometryMagnitudeLimited)                          :: self
    class           (cosmologyFunctionsClass       ), intent(in   ), target   :: cosmologyFunctions_
    double precision                                , intent(in   ), optional :: redshiftMinimum         , redshiftMaximum         , &
         &                                                                       magnitudeApparentMinimum, magnitudeApparentMaximum
    !![
    <constructorAssign variables="magnitudeApparentMinimum, magnitudeApparentMaximum"/>
    !!]

    self%surveyGeometryFullSky=surveyGeometryFullSky(redshiftMinimum=redshiftMinimum,redshiftMaximum=redshiftMaximum,cosmologyFunctions_=cosmologyFunctions_)
    return
  end function magnitudeLimitedConstructorInternal

  subroutine magnitudeLimitedDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryMagnitudeLimited} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryMagnitudeLimited), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine magnitudeLimitedDestructor

  double precision function magnitudeLimitedDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field) result(distanceMinimum)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error                      , only : Error_Report
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    class           (surveyGeometryMagnitudeLimited), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: mass                    , magnitudeAbsolute, &
         &                                                                       luminosity              , starFormationRate
    integer                                         , intent(in   ), optional :: field
    double precision                                                          :: distanceMinimumMagnitude
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(mass             )) call Error_Report(             '`mass` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! If no absolute magnitude is provided, simply return the maximum distance for the survey.
    if (.not.present(magnitudeAbsolute)) then
       distanceMinimum=self%limitDistanceMinimum
       return
    end if
    ! Compute limiting distances. Note that for the magnitudes we have:
    !  m = M₀.₁ + D(z) - 2.5log₁₀(1+z) - K
    ! where D(z)=25+5log(Dₗ) is the regular distance modulus, Dₗ(z)=(1+z)Dᵪ(z) is the luminosity distance, Dᵪ(z) is the comoving
    ! distance, the -2.5log₁₀(1+z) terms accounts for compression of photon frequencies due to redshifting, and K is the
    ! k-correction. Since we do knot know K for each galaxy we neglect it. As such:
    !  D(z) - 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z)] + 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z) √(1+z)] = m - M₀.₁
    ! This is converted to a comoving distance by calling the relevant cosmological conversion function, with the input difference
    ! of magnitudes explicitly noted to contain the -2.5log₁₀(1+z) K-correction factor.
    distanceMinimumMagnitude=self%cosmologyFunctions_%distanceComovingConvert(                                                          &
         &                                                                    output                   =      distanceTypeComoving    , &
         &                                                                    distanceModulusKCorrected=+self%magnitudeApparentMinimum  &
         &                                                                                              -     magnitudeAbsolute         &
         &                                                                   )
    ! Take the larger of the two distances.
    distanceMinimum=max(                                   &
         &                  self%limitDistanceMinimum    , &
         &              min(                               &
         &                  self%limitDistanceMaximum    , &
         &                       distanceMinimumMagnitude  &
         &                 )                               &
         &             )
    return
  end function magnitudeLimitedDistanceMinimum

  double precision function magnitudeLimitedDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field) result(distanceMaximum)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryMagnitudeLimited), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: mass                    , magnitudeAbsolute, &
         &                                                                       luminosity              , starFormationRate
    integer                                         , intent(in   ), optional :: field
    double precision                                                          :: distanceMaximumMagnitude
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(mass             )) call Error_Report(             '`mass` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! If no absolute magnitude is provided, simply return the maximum distance for the survey.
    if (.not.present(magnitudeAbsolute)) then
       distanceMaximum=self%limitDistanceMaximum
       return
    end if
    ! Compute limiting distances. Note that for the magnitudes we have:
    !  m = M₀.₁ + D(z) - 2.5log₁₀(1+z) - K
    ! where D(z)=25+5log(Dₗ) is the regular distance modulus, Dₗ(z)=(1+z)Dᵪ(z) is the luminosity distance, Dᵪ(z) is the comoving
    ! distance, the -2.5log₁₀(1+z) terms accounts for compression of photon frequencies due to redshifting, and K is the
    ! k-correction. Since we do knot know K for each galaxy we neglect it. As such:
    !  D(z) - 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z)] + 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z) √(1+z)] = m - M₀.₁
    ! This is converted to a comoving distance by calling the relevant cosmological conversion function, with the input difference
    ! of magnitudes explicitly noted to contain the -2.5log₁₀(1+z) K-correction factor.
    distanceMaximumMagnitude=self%cosmologyFunctions_%distanceComovingConvert(                                                          &
         &                                                                    output                   =      distanceTypeComoving    , &
         &                                                                    distanceModulusKCorrected=+self%magnitudeApparentMaximum  &
         &                                                                                              -     magnitudeAbsolute         &
         &                                                                   )
    ! Take the smaller of the two distances.
    distanceMaximum=max(                                   &
         &                  self%limitDistanceMinimum    , &
         &              min(                               &
         &                  self%limitDistanceMaximum    , &
         &                       distanceMaximumMagnitude  &
         &             )
    return
  end function magnitudeLimitedDistanceMaximum
