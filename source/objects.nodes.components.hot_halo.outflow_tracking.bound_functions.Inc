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
Contains custom functions for the outflow tracking hot halo component.
!!}

subroutine Node_Component_Hot_Halo_Outflow_Tracking_Mass_Removal_Rate(self,massRemovalRate)
  !!{
  Handle instances where mass is removed from the the outflow tracking hot halo component
  class.
  !!}
  implicit none
  class           (nodeComponentHotHaloOutflowTracking), intent(inout) :: self
  double precision                                     , intent(in   ) :: massRemovalRate

  call self%trackedOutflowMassRate      (-self%trackedOutflowMass      ()*massRemovalRate/self%mass())
  call self%trackedOutflowAbundancesRate(-self%trackedOutflowAbundances()*massRemovalRate/self%mass())
  return
end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Mass_Removal_Rate
