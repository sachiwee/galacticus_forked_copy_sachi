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
  Implements a property extractor class which converts from luminosity to absolute magnitude.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMagnitudeAbsolute">
   <description>A property extractor class which converts from luminosity to absolute magnitude.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMagnitudeAbsolute
     !!{
     A property extractor class which converts from luminosity to absolute magnitude.
     !!}
     private
     class           (nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
     double precision                                       :: offset
     type            (varying_string             )          :: name_                           , description_
   contains
     final     ::                magnitudeAbsoluteDestructor
     procedure :: extract     => magnitudeAbsoluteExtract
     procedure :: type        => magnitudeAbsoluteType
     procedure :: quantity    => magnitudeAbsoluteQuantity
     procedure :: name        => magnitudeAbsoluteName
     procedure :: description => magnitudeAbsoluteDescription
     procedure :: unitsInSI   => magnitudeAbsoluteUnitsInSI
  end type nodePropertyExtractorMagnitudeAbsolute

  interface nodePropertyExtractorMagnitudeAbsolute
     !!{
     Constructors for the \refClass{nodePropertyExtractorMagnitudeAbsolute} output analysis class.
     !!}
     module procedure magnitudeAbsoluteConstructorParameters
     module procedure magnitudeAbsoluteConstructorInternal
  end interface nodePropertyExtractorMagnitudeAbsolute

contains

  function magnitudeAbsoluteConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMagnitudeAbsolute} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorMagnitudeAbsolute)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(nodePropertyExtractorClass            ), pointer       :: nodePropertyExtractor_

    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=nodePropertyExtractorMagnitudeAbsolute(nodePropertyExtractor_)
    class default
       call Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function magnitudeAbsoluteConstructorParameters

  function magnitudeAbsoluteConstructorInternal(nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMagnitudeAbsolute} output analysis property extractor class.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    use :: String_Handling                 , only : String_Upper_Case_First
    implicit none
    type (nodePropertyExtractorMagnitudeAbsolute)                        :: self
    class(nodePropertyExtractorScalar           ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="*nodePropertyExtractor_"/>
    !!]

    ! Find the offset needed to correct for the units of the provided luminosity to the AB system.
    self%offset      =-2.5d0*log10(self%nodePropertyExtractor_%unitsInSI()/luminosityZeroPointAB)
    ! Construct a suitable name and description for the property.
    self%name_       =var_str('magnitudeAbsolute')//String_Upper_Case_First(char(self%nodePropertyExtractor_%name()))
    self%description_=self%nodePropertyExtractor_%description()//' (converted to absolute magnitude)'
    return
  end function magnitudeAbsoluteConstructorInternal

  subroutine magnitudeAbsoluteDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorMagnitudeAbsolute} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMagnitudeAbsolute), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine magnitudeAbsoluteDestructor

  double precision function magnitudeAbsoluteExtract(self,node,instance) result(magnitudeAbsolute)
    !!{
    Implement conversion of luminosity to absolute magnitude.
    !!}
    implicit none
    class           (nodePropertyExtractorMagnitudeAbsolute), intent(inout), target   :: self
    type            (treeNode                              ), intent(inout), target   :: node
    type            (multiCounter                          ), intent(inout), optional :: instance
    double precision                                                                  :: luminosity
    !$GLC attributes unused :: instance

    luminosity=self%nodePropertyExtractor_%extract(node,instance)
    if (luminosity > 0.0d0) then
       magnitudeAbsolute=-2.5d0*log10(luminosity)+self%offset
    else
       ! Zero luminosity, so set infinite absolute magnitude.
       magnitudeAbsolute=huge(0.0d0)
    end if
    return
  end function magnitudeAbsoluteExtract

  function magnitudeAbsoluteType(self) result(type)
    !!{
    Return the type of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeMagnitude
    implicit none
    type (enumerationOutputAnalysisPropertyTypeType)                :: type
    class(nodePropertyExtractorMagnitudeAbsolute   ), intent(inout) :: self
    !$GLC attributes unused :: self

    type=outputAnalysisPropertyTypeMagnitude
    return
  end function magnitudeAbsoluteType

  function magnitudeAbsoluteQuantity(self) result(quantity)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: quantity
    class(nodePropertyExtractorMagnitudeAbsolute       ), intent(inout) :: self
    !$GLC attributes unused :: self

    quantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function magnitudeAbsoluteQuantity

  function magnitudeAbsoluteName(self) result(name)
    !!{
    Return the name of the magnitudeAbsolute property.
    !!}
    implicit none
    type (varying_string                        )                :: name
    class(nodePropertyExtractorMagnitudeAbsolute), intent(inout) :: self

    name=self%name_
    return
  end function magnitudeAbsoluteName

  function magnitudeAbsoluteDescription(self) result(description)
    !!{
    Return a description of the magnitudeAbsolute property.
    !!}
    implicit none
    type (varying_string                        )                :: description
    class(nodePropertyExtractorMagnitudeAbsolute), intent(inout) :: self

    description=self%description_
    return
  end function magnitudeAbsoluteDescription

  double precision function magnitudeAbsoluteUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the absolute magnitude property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorMagnitudeAbsolute), intent(inout) :: self
    !$GLC attributes unused :: self

    unitsInSI=1.0d0
    return
  end function magnitudeAbsoluteUnitsInSI
