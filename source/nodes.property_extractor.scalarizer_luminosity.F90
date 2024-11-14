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

!!{
Contains a module which implements an output analysis property extractor class that scalarizes sum from an array node property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorScalarizerLuminosity">
   <description>An output analysis property extractor class that scalarizes sum of from an array node property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorScalarizerLuminosity
     !!{
     A property extractor output analysis class that scalarizes to return sum of an array node property extractor.
     !!}
     private
     integer                                      :: item                            
     class  (nodePropertyExtractorClass), pointer :: nodePropertyExtractor_ => null()
   contains
     final     ::                scalarizerLuminosityDestructor
     procedure :: extract     => scalarizerLuminosityExtract
     procedure :: quantity    => scalarizerLuminosityQuantity
     procedure :: name        => scalarizerLuminosityName
     procedure :: description => scalarizerLuminosityDescription
     procedure :: unitsInSI   => scalarizerLuminosityUnitsInSI
  end type nodePropertyExtractorScalarizerLuminosity

  interface nodePropertyExtractorScalarizerLuminosity
     !!{
     Constructors for the ``scalarizer'' output analysis class.
     !!}
     module procedure scalarizerLuminosityConstructorParameters
     module procedure scalarizerLuminosityConstructorInternal
  end interface nodePropertyExtractorScalarizerLuminosity

contains

  function scalarizerLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``scalarizer'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodePropertyExtractorScalarizerLuminosity)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (nodePropertyExtractorClass     ), pointer       :: nodePropertyExtractor_
    integer                                                 :: item                  

    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       !![
       <inputParameter>
         <name>item</name>
         <description>The item to scalarize from the array.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    class is (nodePropertyExtractorTuple)
       ! This is as expected.
    class default
       ! "item" is not relevant for non-array extractors.
       item=-1
    end select
    self=nodePropertyExtractorScalarizerLuminosity(item,nodePropertyExtractor_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scalarizerLuminosityConstructorParameters

  function scalarizerLuminosityConstructorInternal(item,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily scalarizer} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (nodePropertyExtractorScalarizerLuminosity)                        :: self
    integer                                 , intent(in   )         :: item                 
    class  (nodePropertyExtractorClass     ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="item, *nodePropertyExtractor_"/>
    !!]
    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       ! This is as expected.
    class is (nodePropertyExtractorTuple)
       ! This is as expected.
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerLuminosityConstructorInternal

  subroutine scalarizerLuminosityDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily scalarizer} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorScalarizerLuminosity), intent(inout) :: self
    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine scalarizerLuminosityDestructor

  double precision function scalarizerLuminosityExtract(self,node,instance)
    !!{
    Implement a scalarizer output analysis.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nodePropertyExtractorScalarizerLuminosity), intent(inout), target         :: self
    type            (treeNode                       ), intent(inout), target         :: node
    type            (multiCounter                   ), intent(inout), optional       :: instance
    class           (nodeComponentBasic             ), pointer                       :: basic
    double precision                                 , allocatable  , dimension(:,:) :: array
    double precision                                 , allocatable  , dimension(  :) :: tuple
   
    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       basic => node%basic()
       if (self%item    > nodePropertyExtractor__%size        (basic%time())) call Error_Report('item exceeds size of array'    //{introspection:location})
       array            =nodePropertyExtractor__%extract(node     ,basic%time   (),instance)
       scalarizerLuminosityExtract=sum(array(self%item, :))
     class is (nodePropertyExtractorTuple)
       basic => node%basic()
       if (self%item    > nodePropertyExtractor__%elementCount(basic%time())) call Error_Report('item exceeds count of tuple'//{introspection:location})
       tuple            =nodePropertyExtractor__%extract(node     ,basic%time   (),instance)
       scalarizerLuminosityExtract=sum(tuple(:))
     class default
       scalarizerLuminosityExtract=0.0d0
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerLuminosityExtract

  function scalarizerLuminosityQuantity(self)
  !!{
    Return the quantity of the scalarizer property.
  !!}
    use :: Error, only : Error_Report
    use :: Output_Analyses_Options, only : enumerationOutputAnalysisPropertyQuantityType
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType  )              :: scalarizerLuminosityQuantity
    class(nodePropertyExtractorScalarizerLuminosity), intent(inout)     :: self
    type (enumerationOutputAnalysisPropertyQuantityType  ) :: quantity

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       scalarizerLuminosityQuantity=nodePropertyExtractor__%quantity()
    class is (nodePropertyExtractorTuple)
       scalarizerLuminosityQuantity=nodePropertyExtractor__%quantity()
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
   end function scalarizerLuminosityQuantity


  function scalarizerLuminosityName(self)
    !!{
    Return the name of the scalarizer property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (varying_string                 )                              :: scalarizerLuminosityName
    class(nodePropertyExtractorScalarizerLuminosity), intent(inout)               :: self
    type (varying_string                 ), allocatable  , dimension(:) :: names
    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       call nodePropertyExtractor__%names(             names)
       scalarizerLuminosityName=names(self%item)
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor__%names(-huge(0.0d0),names)
       scalarizerLuminosityName=names(self%item)
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
   end function scalarizerLuminosityName

  function scalarizerLuminosityDescription(self)
    !!{
    Return a description of the scalarizer property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (varying_string                 )                              :: scalarizerLuminosityDescription
    class(nodePropertyExtractorScalarizerLuminosity), intent(inout)               :: self
    type (varying_string                 ), allocatable  , dimension(:) :: descriptions

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       call nodePropertyExtractor__%descriptions(             descriptions)
       scalarizerLuminosityDescription='Scalarize luminosity emission'
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor__%descriptions(-huge(0.0d0),descriptions)
       scalarizerLuminosityDescription='Scalarize luminosity emission'
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerLuminosityDescription

  double precision function scalarizerLuminosityUnitsInSI(self)
    !!{
    Return the units of the scalarizer property in the SI system.
    !!}
    use :: Error, only : Error_Report
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    class           (nodePropertyExtractorScalarizerLuminosity), intent(inout)               :: self
    double precision                                 , allocatable  , dimension(:) :: unitsInSI

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       unitsInSI          =nodePropertyExtractor__%unitsInSI(            )
       scalarizerLuminosityUnitsInSI=ergs
    class is (nodePropertyExtractorTuple)
       unitsInSI          =nodePropertyExtractor__%unitsInSI(-huge(0.0d0))
       scalarizerLuminosityUnitsInSI=ergs
    class default
       scalarizerLuminosityUnitsInSI=0.0d0
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerLuminosityUnitsInSI
