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
Implements a $1/V_\mathrm{max}$ output analysis weight operator class.
!!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Geometry_Surveys        , only : surveyGeometryClass
  use :: Output_Times            , only : outputTimesClass
  use :: Node_Property_Extractors, only : nodePropertyExtractorClass

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorVolumeMaximum">
   <description>
    An output analysis weight operator class which weights nodes by $1/V_\mathrm{max}$ to mimic what is often done in
    observational analyses. Specifically, the weight of each node is multiplied by $V_j/V_{\mathrm{max}, i}$ where $V_j$ is the
    volume associated with output snapshot $j$, and $V_{\mathrm{max}, i}$ is the maximum volume within which the node could be
    detected. This will up-weight galaxies with small $V_\mathrm{max}$, biasing toward lower redshifts for fainter galaxies, and
    so including any evolution bias that would be present observationally.
   </description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorVolumeMaximum
     !!{
     A cosmological volume corrector analysis weight operator class.
     !!}
     private
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_    => null()
     class           (surveyGeometryClass       ), pointer                   :: surveyGeometry_        => null()
     class           (outputTimesClass          ), pointer                   :: outputTimes_           => null()
     class           (nodePropertyExtractorClass), pointer                   :: nodePropertyExtractor_ => null()
     double precision                            , allocatable, dimension(:) :: volumes
   contains
     final     ::            volumeMaximumDestructor
     procedure :: operate => volumeMaximumOperate
  end type outputAnalysisWeightOperatorVolumeMaximum

  interface outputAnalysisWeightOperatorVolumeMaximum
     !!{
     Constructors for the \refClass{outputAnalysisWeightOperatorVolumeMaximum} output analysis class.
     !!}
     module procedure volumeMaximumConstructorParameters
     module procedure volumeMaximumConstructorInternal
  end interface outputAnalysisWeightOperatorVolumeMaximum

contains

  function volumeMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisWeightOperatorVolumeMaximum} output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputAnalysisWeightOperatorVolumeMaximum)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    class(surveyGeometryClass                      ), pointer       :: surveyGeometry_
    class(outputTimesClass                         ), pointer       :: outputTimes_
    class(nodePropertyExtractorClass               ), pointer       :: nodePropertyExtractor_

    !![
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="surveyGeometry"        name="surveyGeometry_"        source="parameters"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"           source="parameters"/>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisWeightOperatorVolumeMaximum(cosmologyFunctions_,surveyGeometry_,outputTimes_,nodePropertyExtractor_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="surveyGeometry_"       />
    <objectDestructor name="outputTimes_"          />
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function volumeMaximumConstructorParameters

  function volumeMaximumConstructorInternal(cosmologyFunctions_,surveyGeometry_,outputTimes_,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisWeightOperatorVolumeMaximum} output analysis weight operator class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    type            (outputAnalysisWeightOperatorVolumeMaximum)                        :: self
    class           (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class           (surveyGeometryClass                      ), intent(in   ), target :: surveyGeometry_
    class           (outputTimesClass                         ), intent(in   ), target :: outputTimes_
    class           (nodePropertyExtractorClass               ), intent(in   ), target :: nodePropertyExtractor_
    integer         (c_size_t                                 )                        :: iOutput
    double precision                                                                   :: redshiftMinimum       , redshiftMaximum, &
         &                                                                                distanceMinimum       , distanceMaximum
    !![
    <constructorAssign variables="*cosmologyFunctions_, *surveyGeometry_, *outputTimes_, *nodePropertyExtractor_"/>
    !!]

    ! Compute the volume associated with each output.
    allocate(self%volumes(self%outputTimes_%count()))
    do iOutput=1_c_size_t,self%outputTimes_%count()
       write(0,*) "iOutput,1_c_size_t",iOutput==1_c_size_t
       if (iOutput == 1_c_size_t) then
          redshiftMaximum=     self%outputTimes_%redshift(iOutput)
       else
          redshiftMaximum=sqrt(self%outputTimes_%redshift(iOutput)*self%outputTimes_%redshift(iOutput-1_c_size_t))
       end if
       if (iOutput == self%outputTimes_%count()) then
          redshiftMinimum=     self%outputTimes_%redshift(iOutput)
       else
          redshiftMinimum=sqrt(self%outputTimes_%redshift(iOutput)*self%outputTimes_%redshift(iOutput+1_c_size_t))
       end if
       distanceMaximum         =min(                                                                         &
            &                       self%cosmologyFunctions_%distanceComoving           (                    &
            &                       self%cosmologyFunctions_%cosmicTime                  (                   &
            &                       self%cosmologyFunctions_%expansionFactorFromRedshift  (redshiftMaximum)  &
            &                                                                            )                   &
            &                                                                           )                  , &
            &                       self%surveyGeometry_    %distanceMaximum              (               )  &
            &                      )
       distanceMinimum         =max(                                                                         &
            &                       self%cosmologyFunctions_%distanceComoving           (                    &
            &                       self%cosmologyFunctions_%cosmicTime                  (                   &
            &                       self%cosmologyFunctions_%expansionFactorFromRedshift  (redshiftMinimum)  &
            &                                                                            )                   &
            &                                                                           )                  , &
            &                       self%surveyGeometry_    %distanceMinimum              (               )  &
            &                      )
       self%volumes(1_c_size_t)=max(                     &
            &                       +distanceMaximum**3  &
            &                       -distanceMinimum**3, &
            &                       +0.0d0               &
            &                      )
    end do
    if (all(self%volumes <= 0.0d0)) call Error_Report('all output times have zero volume'//{introspection:location})
    return
  end function volumeMaximumConstructorInternal

  subroutine volumeMaximumDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisWeightOperatorVolumeMaximum} output analysis weight operator class.
    !!}
    implicit none
    type(outputAnalysisWeightOperatorVolumeMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%surveyGeometry_"       />
    <objectDestructor name="self%outputTimes_"          />
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine volumeMaximumDestructor

  double precision function volumeMaximumOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex) result(weight)
    !!{
    Implement an volumeMaximum output analysis weight operator.
    !!}
    use            :: Error                  , only : Error_Report
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity, outputAnalysisPropertyQuantityStarFormationRate,outputAnalysisPropertyQuantityMass, outputAnalysisPropertyTypeLinear, &
          &                                           outputAnalysisPropertyTypeMagnitude     , outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisWeightOperatorVolumeMaximum    ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue   , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    double precision                                                               :: distanceMinimum , distanceMaximum       , &
         &                                                                            volumeMaximum
    !$GLC attributes unused :: node

    ! Compute Vₘₐₓ for this node.
    select case (propertyQuantity%ID)
    case (outputAnalysisPropertyQuantityMass      %ID)
       select case (propertyType%ID)
       case (outputAnalysisPropertyTypeLinear   %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   mass             =propertyValue           &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   mass             =propertyValue           &
               &                                                  )
       case (outputAnalysisPropertyTypeLog10    %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   mass             =propertyValueIntrinsic  &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   mass             =propertyValueIntrinsic  &
               &                                                  )
       case default
          distanceMinimum   =+0.0d0
          distanceMaximum   =+0.0d0
          call Error_Report('unsupported property type' //{introspection:location})
       end select
    case (outputAnalysisPropertyQuantityLuminosity%ID)
       select case (propertyType%ID)
       case (outputAnalysisPropertyTypeLinear   %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   luminosity       =propertyValue           &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   luminosity       =propertyValue           &
               &                                                  )
       case (outputAnalysisPropertyTypeLog10    %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   luminosity       =propertyValueIntrinsic  &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   luminosity       =propertyValueIntrinsic  &
               &                                                  )
       case (outputAnalysisPropertyTypeMagnitude    %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   magnitudeAbsolute=propertyValue           &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   magnitudeAbsolute=propertyValue           &
               &                                                  )
       case default
          distanceMinimum   =+0.0d0
          distanceMaximum   =+0.0d0
          call Error_Report('unsupported property type' //{introspection:location})
       end select
    case (outputAnalysisPropertyQuantityStarFormationRate%ID)
       select case (propertyType%ID)
       case (outputAnalysisPropertyTypeLinear   %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   starFormationRate=propertyValue           &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   starFormationRate=propertyValue           &
               &                                                  )
       case (outputAnalysisPropertyTypeLog10    %ID)
          distanceMinimum   =+self%surveyGeometry_%distanceMinimum(                                          &
               &                                                   starFormationRate=propertyValueIntrinsic  &
               &                                                  )
          distanceMaximum   =+self%surveyGeometry_%distanceMaximum(                                          &
               &                                                   starFormationRate=propertyValueIntrinsic  &
               &                                                  )
       case default
          distanceMinimum   =+0.0d0
          distanceMaximum   =+0.0d0
          call Error_Report('unsupported property type' //{introspection:location})
       end select
    case default
       distanceMinimum      =+0.0d0
       distanceMaximum      =+0.0d0
       call    Error_Report('unsupported property class'//{introspection:location})
    end select
    volumeMaximum=+distanceMaximum**3 &
         &        -distanceMinimum**3
    ! Multiply by the volume factor.
    weight=+     weightValue                &
         & *self%volumes      (outputIndex) &
         & /     volumeMaximum
    return
  end function volumeMaximumOperate
