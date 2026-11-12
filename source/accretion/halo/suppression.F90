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

!+    Contributions to this file made by: Sachi Weerasooriya

!!{RST
Contains a module which provides a class that provides a criterion to suppress gas accretion to halos.
!!}
module Accretion_Halo_Suppression
  !!{RST
  Provides a class that provides a criterion to suppress gas accretion to halos.
  !!}
  private

  !![
  <functionClass docformat="rst">
   <name>accretionHaloSuppression</name>
   <descriptiveName>Accretion halo suppression</descriptiveName>
   <description>
   Class providing a criterion to suppress halo accretion based on redshift and virial velocity.
   </description>
   <default>threshold</default>
   <method name="isSuppressed">
    <description>
    Returns true if gas accretion onto the halo from the IGM is suppressed.
    </description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Accretion_Halo_Suppression
