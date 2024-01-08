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
  Implementation of a star formation rate in galactic disks which sets rates in galaxies to zero if time >= quenchTime.
  !!}

  !![
  <starFormationRateDisks name="starFormationRateDisksQuenchTime">
   <description>A star formation rate in galactic disks which sets rates in galaxies to zero if time >= quench time.</description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksQuenchTime
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class(starFormationRateDisksClass), pointer :: starFormationRateDisks_ => null()
     double precision                            :: quenchTimeDisk
   contains
     final     ::         quenchTimeDiskDestructor
     procedure :: rate => quenchTimeSFRateDisks
  end type starFormationRateDisksQuenchTime

  interface starFormationRateDisksQuenchTime
     !!{
     Constructors for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
     !!}
     module procedure quenchTimeDiskConstructorParameters
     module procedure quenchTimeDiskConstructorInternal
  end interface starFormationRateDisksQuenchTime

contains

  function quenchTimeDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationRateDisksQuenchTime)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    double precision                                         :: quenchTimeDisk
    class(starFormationRateDisksClass       ), pointer       :: starFormationRateDisks_

    !![
    <inputParameter>
      <name>quenchTimeDisk</name>
      <defaultValue>14d0</defaultValue>
      <description>The time when SF is QuenchTime (t>=tquench).</description>
      <source>parameters</source>
    </inputParameter>
     !!]

    !![
    <objectBuilder class="starFormationRateDisks" name="starFormationRateDisks_" source="parameters"/>
    !!]
    self=starFormationRateDisksQuenchTime(quenchTimeDisk,starFormationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"/>
    !!]
    return
  end function quenchTimeDiskConstructorParameters

  function quenchTimeDiskConstructorInternal(quenchTimeDisk,starFormationRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type (starFormationRateDisksQuenchTime)                        :: self
    double precision, intent(in   )                                :: quenchTimeDisk
    class(starFormationRateDisksClass       ), intent(in   ), target :: starFormationRateDisks_
    !![
    <constructorAssign variables="quenchTimeDisk, *starFormationRateDisks_"/>
    !!]

    return
  end function quenchTimeDiskConstructorInternal

  subroutine quenchTimeDiskDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateDisksQuenchTime), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_" />
    !!]
    return
  end subroutine quenchTimeDiskDestructor

  double precision function quenchTimeSFRateDisks(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily
    node}. Assumes zero rate for satellites when t>=quench time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(starFormationRateDisksQuenchTime), intent(inout), target   :: self
    type (treeNode                          ), intent(inout), target :: node
    class(nodeComponentBasic                ), pointer               :: basic
    
    basic => node%basic()
    if (basic%time() .ge. self%quenchTimeDisk) then
        quenchTimeSFRateDisks=0.0d0
    else
        quenchTimeSFRateDisks=self%starFormationRateDisks_%rate(node)
       
    end if
    return
  end function quenchTimeSFRateDisks
