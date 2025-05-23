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
Contains a module which provides a class that implements general tasks to be performed by \glc.
!!}

module Tasks
  !!{
  Provides a class that implements general tasks to be performed by \glc.
  !!}
  private

  !![
  <functionClass>
   <name>task</name>
   <descriptiveName>Tasks</descriptiveName>
   <description>Class providing general tasks to be performed by \glc.</description>
   <default>evolveForests</default>
   <functionClassDestroy>no</functionClassDestroy>
   <method name="perform" >
    <description>Perform the task.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>integer, intent(  out), optional :: status</argument>
   </method>
   <method name="requiresOutputFile" >
    <description>Should return true if the task requires the main output file to be open.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     taskRequiresOutputFile=.true.
    </code>
   </method>
  </functionClass>
  !!]

end module Tasks
