!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !+ Contributions to this file made by: Sachi Weerasooriya
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale                    , darkMatterHaloScaleClass
  !![
    <accretionHaloSuppression name="accretionHaloSuppressionThreshold">
    <description>Suppresses accretion when the halo virial velocity is below the quench velocity and the time is later than the quench time.</description>
    </accretionHaloSuppression>
  !!]
  type, extends(accretionHaloSuppressionClass) :: accretionHaloSuppressionThreshold
      private
      class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
      class           (darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_      => null()
      double precision                                                   :: redshiftQuench, timeQuench, velocitySuppressionQuench
    contains
      procedure :: isSuppressed => thresholdIsSuppressed
  end type accretionHaloSuppressionThreshold
  interface accretionHaloSuppressionThreshold
     !!{RST
        Constructors for the :galacticus-class:`accretionHaloSuppressionThreshold` function class.
     !!}
     module procedure accretionHaloSuppressionThresholdConstructorParameters
     module procedure accretionHaloSuppressionThresholdConstructorInternal
  end interface accretionHaloSuppressionThreshold

  contains
    function accretionHaloSuppressionThresholdConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`accretionHaloSuppressionThreshold` function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    implicit none
    type            (accretionHaloSuppressionThreshold)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    double precision                                               :: redshiftQuench, velocitySuppressionQuench, timeQuench
    !![
    <inputParameter>
      <name>redshiftQuench</name>
      <defaultValue>5.0d0</defaultValue>
      <description>Redshift when gas accretion is prevented when redshift is lower than the quench redshift.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocitySuppressionQuench</name>
      <defaultValue>180.0d0</defaultValue>
      <description> scale velocity when gas accretion is prevented.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    timeQuench=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftQuench))
    self=accretionHaloSuppressionThreshold(timeQuench, velocitySuppressionQuench,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    end function accretionHaloSuppressionThresholdConstructorParameters

  function accretionHaloSuppressionThresholdConstructorInternal(timeQuench, velocitySuppressionQuench,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`accretionHaloSuppressionThreshold` gas accretion suppression class.
    !!}
    implicit none
    type            (accretionHaloSuppressionThreshold)                  :: self
    double precision                                     , intent(in   ) :: timeQuench, velocitySuppressionQuench
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="timeQuench, velocitySuppressionQuench,*cosmologyFunctions_,*darkMatterHaloScale_"/>
    !!]
    return
  end function accretionHaloSuppressionThresholdConstructorInternal

  logical function thresholdIsSuppressed(self,node)
    !!{
    Returns true if accretion onto halo is suppressed; time is later than the quench time and halo virial velocity is below the quench velocity.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic,treeNode
    implicit none
    class(accretionHaloSuppressionThreshold), intent(inout) :: self
    type (treeNode                          ), intent(inout):: node
    class(nodeComponentBasic                ), pointer               :: basic
    double precision                                                 :: velocityScale

    basic => node%basic()
    velocityScale=self%darkMatterHaloScale_%velocityVirial(node)
    if (                                                                &
         &  basic%time         (    ) > self%timeQuench                      &
         &   .and.                                                           &
         &  velocityScale < self%velocitySuppressionQuench       &
         & ) then
       thresholdIsSuppressed=.true.
    else
       thresholdIsSuppressed=.false.
    end if
    return
  end function thresholdIsSuppressed


