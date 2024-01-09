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

  !+ Contributions to this file made by: Sachi Weerasooriya
  
  !!{
  Implementation of a star formation rate in galactic disks which sets rates in galaxies to zero at times exceeding a fixed quenching time.
  !!}

  !![
  <starFormationRateDisks name="starFormationRateDisksQuenchTime">
    <description>
      A star formation rate in galactic disks which sets rates in galaxies to zero at times exceeding a fixed quenching
      time. Specifically
      \begin{equation}
      \dot{\phi}_\star(t) \rightarrow \left\{ \begin{array}{ll} \dot{\phi}_\star(t) &amp; \hbox{ if } t &lt; t_\mathrm{quench} \\ 0 &amp; \hbox{ otherwise,} \end{array} \right.
      \end{equation}
      where $t_\mathrm{quench}=${\normalfont \ttfamily [timeQuench]}.
    </description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksQuenchTime
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class           (starFormationRateDisksClass), pointer :: starFormationRateDisks_ => null()
     double precision                                       :: timeQuench
   contains
     final     ::         quenchTimeDestructor
     procedure :: rate => quenchTimeRate
  end type starFormationRateDisksQuenchTime

  interface starFormationRateDisksQuenchTime
     !!{
     Constructors for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
     !!}
     module procedure quenchTimeConstructorParameters
     module procedure quenchTimeConstructorInternal
  end interface starFormationRateDisksQuenchTime

contains

  function quenchTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateDisksQuenchTime)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: timeQuench
    class           (starFormationRateDisksClass     ), pointer       :: starFormationRateDisks_

    !![
    <inputParameter>
      <name>timeQuench</name>
      <defaultValue>14.0d0</defaultValue>
      <description>The time at which star formation should be quenched (truncated to zero).</description>
      <source>parameters</source>
    </inputParameter>
     !!]

    !![
    <objectBuilder class="starFormationRateDisks" name="starFormationRateDisks_" source="parameters"/>
    !!]
    self=starFormationRateDisksQuenchTime(timeQuench,starFormationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"/>
    !!]
    return
  end function quenchTimeConstructorParameters

  function quenchTimeConstructorInternal(timeQuench,starFormationRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type            (starFormationRateDisksQuenchTime)                          :: self
    double precision                                    , intent(in   )         :: timeQuench
    class           (starFormationRateDisksClass       ), intent(in   ), target :: starFormationRateDisks_
    !![
    <constructorAssign variables="timeQuench, *starFormationRateDisks_"/>
    !!]

    return
  end function quenchTimeConstructorInternal

  subroutine quenchTimeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily quenchTime} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateDisksQuenchTime), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_" />
    !!]
    return
  end subroutine quenchTimeDestructor

  double precision function quenchTimeRate(self,node) result(rateStarFormation)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily
    node}. Assumes zero star formation rate at times exceeding the quenching time.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(starFormationRateDisksQuenchTime), intent(inout), target :: self
    type (treeNode                        ), intent(inout), target :: node
    class(nodeComponentBasic              ), pointer               :: basic
    
    basic => node%basic()
    if (basic%time() >= self%timeQuench) then
       rateStarFormation=0.0d0
    else
       rateStarFormation=self%starFormationRateDisks_%rate(node)       
    end if
    return
  end function quenchTimeRate
