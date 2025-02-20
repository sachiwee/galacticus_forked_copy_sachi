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
Contains a module which implements regular expressions by wrapping the GNU C Library implementations.
!!}

! Specify an explicit dependence on the regular_expressions.o object file.
!: $(BUILDPATH)/regular_expressions.o

module Regular_Expressions
  !!{
  Implements regular expressions by wrapping the GNU C Library implementations.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_char, c_int, c_ptr, c_null_ptr
  implicit none
  private
  public :: regEx

  type :: regEx
     !!{
     A regular expression object.
     !!}
     type(c_ptr) :: r=C_NULL_PTR
   contains
     !![
     <methods>
       <method description="Return true if a regular expression matches the supplied {\normalfont \ttfamily string}." method="matches" />
       <method description="Destroy the regex." method="destroy" />
     </methods>
     !!]
     final     ::            Regular_Expression_Destructor
     procedure :: destroy => Regular_Expression_Destroy
     procedure :: matches => Regular_Expression_Match
  end type regEx

  interface regEx
     !!{
     Constructor for regular expression object.
     !!}
     module procedure Regular_Expression_Constructor
  end interface regEx

  interface
     function Regular_Expression_Construct_C(pattern) bind(c,name='Regular_Expression_Construct_C')
       !!{
       Template for a C function that initializes a regular expression.
       !!}
       import
       type     (c_ptr )        :: Regular_Expression_Construct_C
       character(c_char), value :: pattern
     end function Regular_Expression_Construct_C
  end interface

  interface
     subroutine Regular_Expression_Destruct_C(r) bind(c,name='Regular_Expression_Destruct_C')
       !!{
       Template for a C function that destroys a regular expression.
       !!}
       import
       type(c_ptr), value :: r
     end subroutine Regular_Expression_Destruct_C
  end interface

  interface
     function Regular_Expression_Match_C(r,string) bind(c,name='Regular_Expression_Match_C')
       !!{
       Template for a C function that checks for a match with a regular expression.
       !!}
       import
       integer  (c_int )        :: Regular_Expression_Match_C
       type     (c_ptr ), value :: r
       character(c_char), value :: string
     end function Regular_Expression_Match_C
  end interface

contains

  function Regular_Expression_Constructor(regularExpression) result(self)
    !!{
    Constructor for {\normalfont \ttfamily regEx} objects.
    !!}
    implicit none
    type     (regEx)                :: self
    character(len=*), intent(in   ) :: regularExpression

    self%r=Regular_Expression_Construct_C(trim(regularExpression)//char(0))
    return
  end function Regular_Expression_Constructor

  subroutine Regular_Expression_Destructor(self)
    !!{
    Destructor for {\normalfont \ttfamily regEx} objects.
    !!}
    use, intrinsic :: ISO_C_Binding, only : C_Associated
    implicit none
    type(regEx), intent(inout) :: self

    if (C_Associated(self%r)) call Regular_Expression_Destruct_C(self%r)
    self%r=C_NULL_PTR
    return
  end subroutine Regular_Expression_Destructor

  subroutine Regular_Expression_Destroy(self)
    !!{
    Destroy a {\normalfont \ttfamily regEx} object.
    !!}
    implicit none
    class(regEx), intent(inout) :: self

    select type (self)
    type is (regEx)
       call Regular_Expression_Destructor(self)
    end select
    return
  end subroutine Regular_Expression_Destroy

  logical function Regular_Expression_Match(self,string)
    !!{
    Returns true if a {\normalfont \ttfamily regEx} object matches the supplied {\normalfont \ttfamily string}.
    !!}
    implicit none
    class    (regEx)                :: self
    character(len=*), intent(in   ) :: string
    integer                         :: status

    status=Regular_Expression_Match_C(self%r,trim(string)//char(0))
    Regular_Expression_Match=(status==0)
    return
  end function Regular_Expression_Match

end module Regular_Expressions
