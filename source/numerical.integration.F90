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
Contains a module which performs numerical integration.
!!}

! Add dependency on GSL library.
!; gsl

module Numerical_Integration
  !!{
  Implements numerical integration.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_ptr      , c_size_t, c_int, c_double
  use            :: Interface_GSL, only : gslFunction
  implicit none
  private
  public :: integrator

  ! Integrator types.
  !![
  <constant variable="GSL_Integ_Gauss15" gslSymbol="GSL_INTEG_GAUSS15" gslHeader="gsl_integration" type="integer" description="Indicator for 15-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  <constant variable="GSL_Integ_Gauss21" gslSymbol="GSL_INTEG_GAUSS21" gslHeader="gsl_integration" type="integer" description="Indicator for 21-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  <constant variable="GSL_Integ_Gauss31" gslSymbol="GSL_INTEG_GAUSS31" gslHeader="gsl_integration" type="integer" description="Indicator for 31-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  <constant variable="GSL_Integ_Gauss41" gslSymbol="GSL_INTEG_GAUSS41" gslHeader="gsl_integration" type="integer" description="Indicator for 41-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  <constant variable="GSL_Integ_Gauss51" gslSymbol="GSL_INTEG_GAUSS51" gslHeader="gsl_integration" type="integer" description="Indicator for 51-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  <constant variable="GSL_Integ_Gauss61" gslSymbol="GSL_INTEG_GAUSS61" gslHeader="gsl_integration" type="integer" description="Indicator for 61-point Gauss-Kronrod integration rule." reference="GSL" referenceURL="https://www.gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_qag" group="GSL"/>
  !!]

  !![
  <stateStorable class="integrator">
   <integrator>
    <methodCall method="GSLReallocate" arguments="gslFree=.true." />
   </integrator>
  </stateStorable>
  !!]

  !![
  <deepCopyActions class="integrator">
    <integrator>
      <methodCall method="GSLReallocate" arguments="gslFree=.false."/>
    </integrator>
  </deepCopyActions>
  !!]
  
  type :: integrator
     !!{
     Class for performing numerical integrations.
     !!}
     private
     type            (c_ptr            )        , allocatable :: integrandFunction, integrationWorkspace
     procedure       (integrandTemplate), nopass, pointer     :: integrand
     integer                                                  :: integrationRule
     integer         (c_size_t         )                      :: intervalsMaximum
     double precision                                         :: toleranceAbsolute, toleranceRelative
     logical                                                  :: hasSingularities
   contains
     !![
     <methods>
       <method description="Evaluate the integral."                    method="integrate"    />
       <method description="Set tolerances to use in this integrator." method="toleranceSet" />
       <method description="Allocate GSL objects."                     method="gslAllocate"  />
       <method description="Reallocate GSL objects."                   method="gslReallocate"/>
     </methods>
     !!]
     final     ::                  integratorDestructor
     procedure :: integrate     => integratorIntegrate
     procedure :: toleranceSet  => integratorToleranceSet
     procedure :: gslAllocate   => integratorGSLAllocate
     procedure :: gslReallocate => integratorGSLReallocate
  end type integrator

  interface integrator
     !!{
     Interface to constructor for integrators.
     !!}
     module procedure :: integratorConstructor
  end interface integrator
  
  interface
     !!{
     Interfaces to GSL integration functions.
     !!}
     function gsl_integration_qag(f,a,b,epsabs,epsrel,limit,key,workspace,result,abserr) bind(c,name='gsl_integration_qag')
       !!{
       Template for the GSL QAG integration function.
       !!}
       import c_ptr, c_size_t, c_int, c_double
       integer(c_int   )                       :: gsl_integration_qag
       type   (c_ptr   )               , value :: f
       real   (c_double), intent(in   ), value :: a                  , b     , &
            &                                     epsabs             , epsrel
       integer(c_size_t), intent(in   ), value :: limit
       integer(c_int   ), intent(in   ), value :: key
       type   (c_ptr   ), intent(in   ), value :: workspace
       real   (c_double), intent(  out)        :: result             , abserr
     end function gsl_integration_qag

     function gsl_integration_qags(f,a,b,epsabs,epsrel,limit,workspace,result,abserr) bind(c,name='gsl_integration_qags')
       !!{
       Template for the GSL QAGS integration function.
       !!}
       import c_ptr, c_size_t, c_int, c_double
       integer(c_int   )                       :: gsl_integration_qags
       type   (c_ptr   )               , value :: f
       real   (c_double), intent(in   ), value :: a                   , b     , &
            &                                     epsabs              , epsrel
       integer(c_size_t), intent(in   ), value :: limit
       type   (c_ptr   ), intent(in   ), value :: workspace
       real   (c_double), intent(  out)        :: result              , abserr
     end function gsl_integration_qags

     function gsl_integration_workspace_alloc(n) bind(c,name='gsl_integration_workspace_alloc')
       !!{
       Template for GSL integration workspace allocation function.
       !!}
       import c_ptr, c_size_t
       type   (c_ptr   )                       :: gsl_integration_workspace_alloc
       integer(c_size_t), intent(in   ), value :: n
     end function gsl_integration_workspace_alloc

     subroutine gsl_integration_workspace_free(w) bind(c,name='gsl_integration_workspace_free')
       !!{
       Template for GSL integration workspace deallocation function.
       !!}
       import c_ptr
       type(c_ptr), intent(in   ), value :: w
     end subroutine gsl_integration_workspace_free
  end interface

  ! Integrand interface.
  abstract interface
     double precision function integrandTemplate(x)
       double precision, intent(in   ) :: x
     end function integrandTemplate
  end interface

  ! Integrand function.
  procedure(integrandTemplate), pointer :: currentIntegrand
  !$omp threadprivate(currentIntegrand)

contains

  function integratorConstructor(integrand,toleranceAbsolute,toleranceRelative,intervalsMaximum,hasSingularities,integrationRule) result(self)
    !!{
    Constructor for {\normalfont \ttfamily integrator} objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (integrator       )                          :: self
    procedure       (integrandTemplate)                          :: integrand
    integer                            , intent(in   ), optional :: integrationRule
    integer         (c_size_t         ), intent(in   ), optional :: intervalsMaximum
    double precision                   , intent(in   ), optional :: toleranceAbsolute, toleranceRelative
    logical                            , intent(in   ), optional :: hasSingularities
    !![
    <optionalArgument name="integrationRule"   defaultsTo="GSL_Integ_Gauss61"/>
    <optionalArgument name="intervalsMaximum"  defaultsTo="1000_c_size_t"    />
    <optionalArgument name="hasSingularities"  defaultsTo=".false."          />
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0"            />
    <optionalArgument name="toleranceRelative" defaultsTo="0.0d0"            />
    !!]

    ! Validate input.
    if     (                             &
         &   toleranceAbsolute_ <= 0.0d0 &
         &  .and.                        &
         &   toleranceRelative_ <= 0.0d0 &
         & ) call Error_Report('at least one of absolute or relative tolerance must be greater than zero'//{introspection:location})
    self%integrand         => integrand
    self%intervalsMaximum  =  intervalsMaximum_
    self%hasSingularities  =  hasSingularities_
    self%integrationRule   =  integrationRule_
    self%toleranceAbsolute =  toleranceAbsolute_
    self%toleranceRelative =  toleranceRelative_
    call self%gslAllocate()
    return
  end function integratorConstructor
  
  subroutine integratorDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily integrator} objects.
    !!}
    use :: Interface_GSL, only : gslFunctionDestroy
    implicit none
    type(integrator), intent(inout) :: self

    if (allocated(self%integrandFunction   )) then
       call gslFunctionDestroy            (self%integrandFunction   )
       deallocate(self%integrandFunction   )
    end if
    if (allocated(self%integrationWorkspace)) then
       call gsl_integration_workspace_free(self%integrationWorkspace)
       deallocate(self%integrationWorkspace)
    end if
    return
  end subroutine integratorDestructor
  
  subroutine integratorGSLAllocate(self)
    !!{
    Allocate GSL objects.
    !!}
    implicit none
    class(integrator), intent(inout) :: self

    if (.not.allocated(self%integrationWorkspace)) then
       allocate(self%integrationWorkspace)
       self%integrationWorkspace=gsl_integration_workspace_alloc(self%intervalsMaximum)
    end if
    if (.not.allocated(self%integrandFunction)) then
       allocate(self%integrandFunction)
       self%integrandFunction   =gslFunction                    (     integrandWrapper)
    end if
    return
  end subroutine integratorGSLAllocate

  subroutine integratorGSLReallocate(self,gslFree)
    !!{
    Reallocate GSL objects.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_null_ptr
    use            :: Interface_GSL, only : gslFunctionDestroy
    implicit none
    class  (integrator), intent(inout) :: self
    logical            , intent(in   ) :: gslFree
    
    if (allocated(self%integrationWorkspace)) then
       if (gslFree) then
          call gsl_integration_workspace_free(self%integrationWorkspace)
       else
          self%integrationWorkspace=c_null_ptr
       end if
       deallocate(self%integrationWorkspace)
    end if
    if (allocated(self%integrandFunction   )) then
       if (gslFree) then
          call gslFunctionDestroy            (self%integrandFunction   )
       else
          self%integrandFunction   =c_null_ptr
       end if
       deallocate(self%integrandFunction   )
    end if
    call self%GSLAllocate()
    return
  end subroutine integratorGSLReallocate
  
  subroutine integratorToleranceSet(self,toleranceAbsolute,toleranceRelative)
    !!{
    Reset tolerance for {\normalfont \ttfamily integrator} objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (integrator)                          :: self
    double precision            , intent(in   ), optional :: toleranceAbsolute, toleranceRelative
    !![
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0"/>
    <optionalArgument name="toleranceRelative" defaultsTo="0.0d0"/>
    !!]

    ! Validate input.
    if     (                             &
         &   toleranceAbsolute_ <= 0.0d0 &
         &  .and.                        &
         &   toleranceRelative_ <= 0.0d0 &
         & ) call Error_Report('at least one of absolute or relative tolerance must be greater than zero'//{introspection:location})
    self%toleranceAbsolute=toleranceAbsolute_
    self%toleranceRelative=toleranceRelative_    
    return
  end subroutine integratorToleranceSet
  
  recursive double precision function integratorIntegrate(self,limitLower,limitUpper,status)
    !!{
    Perform a numerical integration.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funptr
    use            :: Error        , only : errorStatusSuccess, GSL_Error_Handler_Abort_Off, GSL_Error_Handler_Abort_On
    implicit none
    class           (integrator       ), intent(inout)           :: self
    double precision                   , intent(in   )           :: limitLower       , limitUpper
    integer                            , intent(  out), optional :: status
    procedure       (integrandTemplate), pointer                 :: previousIntegrand
    integer                                                      :: status_
    double precision                                             :: errorAbsolute

    ! Store a pointer to the current integrand (so that we can restore it later), and set the current integrand to our integrand.
    previousIntegrand => currentIntegrand
    currentIntegrand  => self            %integrand
    ! Set error handler if necessary.
    if (present(status)) then
       call GSL_Error_Handler_Abort_Off()
       status_=errorStatusSuccess
    end if
    ! Do the integration
    if (self%hasSingularities) then
       status_=gsl_integration_qags(                           &
            &                       self%integrandFunction   , &
            &                            limitLower          , &
            &                            limitUpper          , &
            &                       self%toleranceAbsolute   , &
            &                       self%toleranceRelative   , &
            &                       self%intervalsMaximum    , &
            &                       self%integrationWorkspace, &
            &                            integratorIntegrate , &
            &                            errorAbsolute         &
            &                      )
    else
       status_=gsl_integration_qag (                           &
            &                       self%integrandFunction   , &
            &                            limitLower          , &
            &                            limitUpper          , &
            &                       self%toleranceAbsolute   , &
            &                       self%toleranceRelative   , &
            &                       self%intervalsMaximum    , &
            &                       self%integrationRule     , &
            &                       self%integrationWorkspace, &
            &                            integratorIntegrate , &
            &                            errorAbsolute         &
            &                      )
    end if
    ! Reset error handler.
    if (present(status)) then
       status=status_
       call GSL_Error_Handler_Abort_On()
    end if
    ! Restore the previous integrand.
    currentIntegrand => previousIntegrand
    return
  end function integratorIntegrate

  function integrandWrapper(x) bind(c)
    !!{
    Wrapper function used for \gls{gsl} integration functions.
    !!}
    implicit none
    real(c_double)                       :: integrandWrapper
    real(c_double), intent(in   ), value :: x

    integrandWrapper=currentIntegrand(x)
    return
  end function integrandWrapper

end module Numerical_Integration
