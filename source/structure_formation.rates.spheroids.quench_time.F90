!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of a star formation rate in galactic disks which sets rates in
galaxies to zero if time >= quenchTimeSpheroid.
  !!}

  !![
  <starFormationRateSpheroids name="starFormationRateSpheroidsQuenchTime">
   <description>A star formation rate in galactic disks which sets rates in galaxies to zero if time >= quench time.</description>
  </starFormationRateSpheroids>
  !!]
  type, extends(starFormationRateSpheroidsClass) :: starFormationRateSpheroidsQuenchTime
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class(starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     double precision                            :: quenchTimeSpheroid
   contains
     final     ::         quenchTimeSpheroidDestructor
     procedure :: rate => quenchTimeSFRateSpheroids
  end type starFormationRateSpheroidsQuenchTime

  interface starFormationRateSpheroidsQuenchTime
     !!{
     Constructors for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
     !!}
     module procedure quenchTimeSpheroidConstructorParameters
     module procedure quenchTimeSpheroidConstructorInternal
  end interface starFormationRateSpheroidsQuenchTime

contains

  function quenchTimeSpheroidConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationRateSpheroidsQuenchTime)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    double precision                                         :: quenchTimeSpheroid
    class(starFormationRateSpheroidsClass       ), pointer       :: starFormationRateSpheroids_

    !![
    <inputParameter>
      <name>quenchTimeSpheroid</name>
      <defaultValue>14d0</defaultValue>
      <description>The time when SF is QuenchTime (t>=tquench).</description>
      <source>parameters</source>
    </inputParameter>
     !!]

    !![
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=starFormationRateSpheroidsQuenchTime(quenchTimeSpheroid,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function quenchTimeSpheroidConstructorParameters

  function quenchTimeSpheroidConstructorInternal(quenchTimeSpheroid,starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type (starFormationRateSpheroidsQuenchTime)                        :: self
    double precision, intent(in   )                                :: quenchTimeSpheroid
    class(starFormationRateSpheroidsClass       ), intent(in   ), target :: starFormationRateSpheroids_
    !![
    <constructorAssign variables="quenchTimeSpheroid, *starFormationRateSpheroids_"/>
    !!]

    return
  end function quenchTimeSpheroidConstructorInternal

  subroutine quenchTimeSpheroidDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateSpheroidsQuenchTime), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSpheroids_" />
    !!]
    return
  end subroutine quenchTimeSpheroidDestructor

  double precision function quenchTimeSFRateSpheroids(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily
    node}. Assumes zero rate for satellites when t>=quench time spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(starFormationRateSpheroidsQuenchTime), intent(inout), target   :: self
    type (treeNode                          ), intent(inout), target :: node
    class(nodeComponentBasic                ), pointer               :: basic
    
    basic => node%basic()
    if (basic%time() .ge. self%quenchTimeSpheroid) then
        quenchTimeSFRateSpheroids=0.0d0
    else
        quenchTimeSFRateSpheroids=self%starFormationRateSpheroids_%rate(node)
       
    end if
    return
  end function quenchTimeSFRateSpheroids
