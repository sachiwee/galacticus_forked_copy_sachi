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
!+    Contributions by Sachi Weerasooriya
!!{
Implements the geometry of the SDSS survey used by \cite{favole_galaxy_2024}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryFavole2024SDSS">
   <description>Implements the geometry of the SDSS survey of \cite{strauss_galaxy_2022}.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryFavole2024SDSS
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                    favole2024SDSSDestructor
     procedure :: distanceMinimum => favole2024SDSSDistanceMinimum
     procedure :: distanceMaximum => favole2024SDSSDistanceMaximum
  end type surveyGeometryFavole2024SDSS

  interface surveyGeometryFavole2024SDSS
     !!{
     Constructors for the \cite{favole_galaxy_2024} survey geometry class.
     !!}
     module procedure favole2024SDSSConstructorParameters
     module procedure favole2024SDSSConstructorInternal
  end interface surveyGeometryFavole2024SDSS

contains

  function favole2024SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \cite{favole_galaxy_2024} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (surveyGeometryFavole2024SDSS)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=surveyGeometryFavole2024SDSS(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function favole2024SDSSConstructorParameters

  function favole2024SDSSConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Default constructor for the \cite{favole_galaxy_2024} survey geometry class.
    !!}
    implicit none
    type (surveyGeometryFavole2024SDSS)                        :: self
    class(cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    ! Initialize geometry.
    call self%initialize()
    return
  end function favole2024SDSSConstructorInternal

  subroutine favole2024SDSSDestructor(self)
    !!{
    Destructor for the ``favole2024SDSS'' survey geometry class.
    !!}
    implicit none
    type(surveyGeometryFavole2024SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine favole2024SDSSDestructor

  double precision function favole2024SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryFavole2024SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity, starFormationRate
    integer                                             , intent(in   ), optional :: field
    !$GLC attributes unused :: field, mass, luminosity, starFormationRate, magnitudeAbsolute

    ! Compute limiting distances. This is due only to the redshift limit.
   favole2024SDSSDistanceMinimum=self   %cosmologyFunctions_%distanceComoving           (          &
         &                               self  %cosmologyFunctions_%cosmicTime                  (         &
         &                                self %cosmologyFunctions_%expansionFactorFromRedshift  (        &
         &                                                                                        +1.0d-3 &
         &                                                                                       )        &
         &                                                                                      )         &
         &                                                                                     )
    return
  end function favole2024SDSSDistanceMinimum

  double precision function favole2024SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options     , only : distanceTypeComoving
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    class           (surveyGeometryFavole2024SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass                           , magnitudeAbsolute        , &
         &                                                                           luminosity, starFormationRate
    integer                                             , intent(in   ), optional :: field
    double precision                                    , parameter               :: fluxLimiting           =2.0d-16 ! W m⁻².
    double precision                                                              :: distanceMaximumRedshift        , distanceMaximumLuminosity, &
         &                                                                           distanceLuminosity
    !$GLC attributes unused :: field, mass, magnitudeAbsolute, starFormationRate

    ! Validate input.
    
    if (.not.present(luminosity)) call Error_Report('luminosity must be supplied '//{introspection:location})
    ! Compute limiting distances. We find the luminosity distance from the supplied luminosity and the limiting flux of the survey.
    distanceLuminosity      =sqrt(                 &
         &                        +luminosity      &
         &                        *ergs            &
         &                        /4.0d0           &
         &                        /Pi              &
         &                        /fluxLimiting    &
         &                        /megaParsec  **2 &
         &                       )
    distanceMaximumRedshift =self   %cosmologyFunctions_%distanceComoving           (                                         &
         &                    self  %cosmologyFunctions_%cosmicTime                  (                                        &
         &                     self %cosmologyFunctions_%expansionFactorFromRedshift  (                                       &
         &                                                                                              1.0d-1                &
         &                                                                            )                                       &
         &                                                                           )                                        &
         &                                                                          )
    distanceMaximumLuminosity=self   %cosmologyFunctions_%distanceComovingConvert    (                                        &
         &                                                                           output            =distanceTypeComoving, &
         &                                                                           distanceLuminosity=distanceLuminosity    &
         &                                                                          )
    ! Take the smaller of the two distances.
    favole2024SDSSDistanceMaximum=min(                           &
         &                                  distanceMaximumRedshift  , &
         &                                  distanceMaximumLuminosity  &
         &                                 )
    return
  end function favole2024SDSSDistanceMaximum
