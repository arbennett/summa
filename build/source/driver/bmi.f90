module bmi_summa
    
    use nrtype                                                  ! variable types, etc.
    use netcdf                                                  ! netcdf libraries
    ! provide access to subroutines and functions
    use summaFileManager,       only: summa_SetDirsUndPhiles    ! sets directories and filenames
    use module_sf_noahmplsm,    only: read_mp_veg_parameters    ! module to read NOAH vegetation tables
    use module_sf_noahmplsm,    only: redprm                    ! module to assign more Noah-MP parameters
    use module_sf_noahmplsm,    only: isWater                   ! parameter for water land cover type
    use nr_utility_module,      only: arth                      ! get a sequence of numbers
    use ascii_util_module,      only: file_open                 ! open ascii file
    use ascii_util_module,      only: get_vlines                ! read a vector of non-comment lines from an ASCII file
    use ascii_util_module,      only: split_line                ! extract the list of variable names from the character string
    use time_utils_module,      only: elapsedSec                ! calculate the elapsed time
    use allocspace_module,      only: allocGlobal               ! module to allocate space for global data structures
    use allocspace_module,      only: allocLocal                ! module to allocate space for local data structures
    use childStruc_module,      only: childStruc                ! module to create a child data structure
    use mDecisions_module,      only: mDecisions                ! module to read model decisions
    use popMetadat_module,      only: popMetadat                ! module to populate metadata structures
    use flxMapping_module,      only: flxMapping                ! module to map fluxes to states
    use checkStruc_module,      only: checkStruc                ! module to check metadata structures
    use def_output_module,      only: def_output                ! module to define model output
    use ffile_info_module,      only: ffile_info                ! module to read information on forcing datafile
    use read_attrb_module,      only: read_dimension            ! module to read dimensions of GRU and HRU
    use read_attrb_module,      only: read_attrb                ! module to read local attributes
    use read_pinit_module,      only: read_pinit                ! module to read initial model parameter values
    use paramCheck_module,      only: paramCheck                ! module to check consistency of model parameters
    use check_icond_module,     only: check_icond               ! module to check initial conditions
    use read_icond_module,      only: read_icond                ! module to read initial conditions
    use read_icond_module,      only: read_icond_nlayers        ! module to read initial conditions
    use pOverwrite_module,      only: pOverwrite                ! module to overwrite default parameter values with info from the Noah tables
    use read_param_module,      only: read_param                ! module to read model parameter sets
    use ConvE2Temp_module,      only: E2T_lookup                ! module to calculate a look-up table for the temperature-enthalpy conversion
    use var_derive_module,      only: calcHeight                ! module to calculate height at layer interfaces and layer mid-point
    use var_derive_module,      only: v_shortcut                ! module to calculate "short-cut" variables
    use var_derive_module,      only: rootDensty                ! module to calculate the vertical distribution of roots
    use var_derive_module,      only: satHydCond                ! module to calculate the saturated hydraulic conductivity in each soil layer
    use var_derive_module,      only: fracFuture                ! module to calculate the fraction of runoff in future time steps (time delay histogram)
    use read_force_module,      only: read_force                ! module to read model forcing data
    use derivforce_module,      only: derivforce                ! module to compute derived forcing data
    use modelwrite_module,      only: writeParm,writeTime       ! module to write model attributes and parameters
    use modelwrite_module,      only: writeData,writeBasin      ! module to write model output
    use modelwrite_module,      only: writeRestart              ! module to write model Restart
    use vegPhenlgy_module,      only: vegPhenlgy                ! module to compute vegetation phenology
    use coupled_em_module,      only: coupled_em                ! module to run the coupled energy and mass model
    use groundwatr_module,      only: groundwatr                ! module to simulate regional groundwater balance
    use qTimeDelay_module,      only: qOverland                 ! module to route water through an "unresolved" river network
    use netcdf_util_module,     only: nc_file_close             ! module to handle netcdf stuff for inputs and outputs
    ! provide access to file paths
    use summaFileManager,       only: SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
    use summaFileManager,       only: MODEL_INITCOND                    ! name of model initial conditions file
    use summaFileManager,       only: LOCAL_ATTRIBUTES                  ! name of model initial attributes file
    use summaFileManager,       only: OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
    use summaFileManager,       only: LOCALPARAM_INFO,BASINPARAM_INFO   ! files defining the default values and constraints for model parameters
    ! provide access to the derived types to define the data structures
    use data_types, only:&
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (dp)
    use data_types,             only: extended_info          ! extended metadata structure
    ! provide access to runtime options
    use globalData,             only: iRunModeFull, iRunModeGRU, iRunModeHRU
    ! provide access to metadata structures
    use globalData, only: time_meta, forc_meta, attr_meta, type_meta    ! metadata structures
    use globalData, only: prog_meta, diag_meta, flux_meta               ! metadata structures
    use globalData, only: mpar_meta, indx_meta                          ! metadata structures
    use globalData, only: bpar_meta, bvar_meta                          ! metadata structures
    use globalData, only: averageFlux_meta                              ! metadata for time-step average fluxes
    use globalData, only: model_decisions                               ! model decision structure
    ! provide access to global data
    use globalData, only: refTime                                       ! reference time
    use globalData, only: startTime                                     ! start time
    use globalData, only: finshTime                                     ! end time
    use globalData, only: doJacobian                                    ! flag to compute the Jacobian
    use globalData, only: gru_struc                                     ! gru-hru mapping structures
    use globalData, only: localParFallback                              ! local column default parameters
    use globalData, only: basinParFallback                              ! basin-average default parameters
    use globalData, only: structInfo                                    ! information on the data structures                  
    use globalData, only: numtim                                        ! number of time steps
    use globalData, only: urbanVegCategory                              ! vegetation category for urban areas
    use globalData, only: globalPrintFlag                               ! global print flag
    use globalData, only: integerMissing                                ! missing integer value
    ! provide access to Noah-MP parameters
    use NOAHMP_VEG_PARAMETERS,  only: SAIM, LAIM                            ! 2-d tables for stem area index and leaf area index (vegType,month)
    use NOAHMP_VEG_PARAMETERS,  only: HVT, HVB                              ! height at the top and bottom of vegetation (vegType)
    use noahmp_globals,         only: RSMIN                                 ! minimum stomatal resistance (vegType)
    use var_lookup,             only: maxvarForc, maxvarProg, maxvarDiag    ! size of variable vectors 
    use var_lookup,             only: maxvarFlux, maxvarIndx, maxvarBvar    ! size of variable vectors 
    ! provide access to the named variables that describe elements of parent model structures
    use var_lookup, only: iLookTIME, iLookFORCE                         ! look-up values for time and forcing data structures
    use var_lookup, only: iLookTYPE                                     ! look-up values for classification of veg, soils etc.
    use var_lookup, only: iLookATTR                                     ! look-up values for local attributes
    use var_lookup, only: iLookPARAM                                    ! look-up values for local column model parameters
    use var_lookup, only: iLookINDEX                                    ! look-up values for local column index variables
    use var_lookup, only: iLookPROG                                     ! look-up values for local column model prognostic (state) variables
    use var_lookup, only: iLookDIAG                                     ! look-up values for local column model diagnostic variables 
    use var_lookup, only: iLookFLUX                                     ! look-up values for local column model fluxes 
    use var_lookup, only: iLookBVAR                                     ! look-up values for basin-average model variables
    use var_lookup, only: iLookBPAR                                     ! look-up values for basin-average model parameters
    use var_lookup, only: iLookDECISIONS                                ! look-up values for model decisions
    use var_lookup, only: iLookVarType                                  ! look-up values for variable type structure
    ! provide access to the named variables that describe elements of child  model structures
    use var_lookup, only: childFLUX_MEAN                                ! look-up values for timestep-average model fluxes
    ! provide access to the named variables that describe model decisions
    use mDecisions_module,  only:  &                                    ! look-up values for method used to compute derivative
                            numerical,    &                             ! numerical solution
                            analytical,   &                             ! analytical solution
                            monthlyTable, &                             ! LAI/SAI taken directly from a monthly table for different vegetation classes
                            specified,    &                             ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
                            localColumn,  &                             ! separate groundwater representation in each local soil column
                            singleBasin                                 ! single groundwater store over the entire basin
    use output_stats,       only: calcStats                             ! module for compiling output statistics
    use globalData,         only: nFreq, outFreq                        ! model output files
    use globalData,         only: ncid                                  ! file id of netcdf output file
    use var_lookup,         only: maxFreq                               ! maximum # of output files

implicit none

    ! *****************************************************************************
    ! BMI codes
    ! *****************************************************************************

    integer, parameter :: BMI_VAR_TYPE_UNKNOWN = 0
    integer, parameter :: BMI_VAR_TYPE_CHAR = 1
    integer, parameter :: BMI_VAR_TYPE_UNSIGNED_CHAR = 2
    integer, parameter :: BMI_VAR_TYPE_INT = 3
    integer, parameter :: BMI_VAR_TYPE_LONG = 4
    integer, parameter :: BMI_VAR_TYPE_UNSIGNED_INT = 5
    integer, parameter :: BMI_VAR_TYPE_UNSIGNED_LONG = 6
    integer, parameter :: BMI_VAR_TYPE_FLOAT = 7
    integer, parameter :: BMI_VAR_TYPE_DOUBLE = 8
    integer, parameter :: BMI_VAR_TYPE_NUMBER = 9

    integer, parameter :: BMI_GRID_TYPE_UNKNOWN = 0
    integer, parameter :: BMI_GRID_TYPE_UNIFORM = 1
    integer, parameter :: BMI_GRID_TYPE_RECTILINEAR = 2
    integer, parameter :: BMI_GRID_TYPE_STRUCTURED = 3
    integer, parameter :: BMI_GRID_TYPE_UNSTRUCTURED = 4
    integer, parameter :: BMI_GRID_TYPE_NUMBER = 5

    integer, parameter :: BMI_MAXVARNAMESTR = 31
    integer, parameter :: BMI_MAXCOMPNAMESTR = 31
    integer, parameter :: BMI_MAXUNITSSTR = 31

    integer, parameter :: BMI_CHAR = 1
    integer, parameter :: BMI_UNSIGNED_CHAR = 1
    integer, parameter :: BMI_INT = 2
    integer, parameter :: BMI_LONG = 4
    integer, parameter :: BMI_UNSIGNED_INT = 2
    integer, parameter :: BMI_UNSIGNED_LONG = 4
    integer, parameter :: BMI_FLOAT = 4
    integer, parameter :: BMI_DOUBLE = 8

    integer, parameter :: BMI_FAILURE = 0
    integer, parameter :: BMI_SUCCESS = 1

    ! *****************************************************************************
    ! (0) variable definitions
    ! *****************************************************************************

    type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
    type(gru_hru_doubleVec)          :: progStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_doubleVec)          :: diagStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_doubleVec)          :: fluxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
    type(gru_hru_doubleVec)          :: indxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model indices
    type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat        -- basin-average variabl
    ! define the primary data structures (scalars)
    type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
    type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
    type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
    type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
    ! define the primary data structures (variable length vectors)
    type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
    type(gru_hru_doubleVec)          :: mparStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
    type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
    ! define the basin-average structures
    type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
    type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables
    ! define the ancillary data structures
    type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:)     -- default model parameters
    ! define indices
    integer(i4b)                     :: iStruct                    ! loop through data structures
    integer(i4b)                     :: iGRU
    integer(i4b)                     :: iHRU,jHRU,kHRU             ! index of the hydrologic response unit
    integer(i4b)                     :: nGRU                       ! number of grouped response units
    integer(i4b)                     :: nHRU                       ! number of global hydrologic response units
    integer(i4b)                     :: hruCount                   ! number of local hydrologic response units
    integer(i4b)                     :: modelTimeStep=0            ! index of model time step
    integer(i4b)                     :: waterYearTimeStep=0        ! index of water year
    integer(i4b),dimension(maxFreq)  :: outputTimeStep=0           ! timestep in output files
    ! define the time output
    logical(lgt)                     :: printProgress              ! flag to print progress
    integer(i4b),parameter           :: ixProgress_im=1000         ! named variable to print progress once per month
    integer(i4b),parameter           :: ixProgress_id=1001         ! named variable to print progress once per day
    integer(i4b),parameter           :: ixProgress_ih=1002         ! named variable to print progress once per hour
    integer(i4b),parameter           :: ixProgress_never=1003      ! named variable to print progress never
    integer(i4b)                     :: ixProgress=ixProgress_id   ! define frequency to write progress
    ! define the re-start file
    logical(lgt)                     :: printRestart               ! flag to print a re-start file
    integer(i4b),parameter           :: ixRestart_iy=1000          ! named variable to print a re-start file once per year
    integer(i4b),parameter           :: ixRestart_im=1001          ! named variable to print a re-start file once per month
    integer(i4b),parameter           :: ixRestart_id=1002          ! named variable to print a re-start file once per day
    integer(i4b),parameter           :: ixRestart_never=1003       ! named variable to print a re-start file never
    integer(i4b)                     :: ixRestart=ixRestart_never  ! define frequency to write restart files
    ! define output file
    integer(i4b)                     :: ctime1(8)                  ! initial time
    character(len=256)               :: output_fileSuffix=''       ! suffix for the output file
    character(len=256)               :: summaFileManagerFile=''    ! path/name of file defining directories and files
    character(len=256)               :: fileout=''                 ! output filename
    ! define model control structures
    integer(i4b)                     :: nLayers                    ! total number of layers
    integer(i4b),parameter           :: no=0                       ! .false.
    integer(i4b),parameter           :: yes=1                      ! .true.
    logical(lgt)                     :: computeVegFluxFlag         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow) 
    type(hru_i),allocatable          :: computeVegFlux(:)          ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow) 
    type(hru_d),allocatable          :: dt_init(:)                 ! used to initialize the length of the sub-step for each HRU
    type(hru_d),allocatable          :: upArea(:)                  ! area upslope of each HRU 
    ! general local variables        
    integer(i4b)                     :: ivar                       ! index of model variable
    real(dp)                         :: fracHRU                    ! fractional area of a given HRU (-)
    logical(lgt)                     :: flux_mask(maxvarFlux)      ! mask defining desired flux variables
    integer(i4b)                     :: forcNcid=integerMissing    ! netcdf id for current netcdf forcing file
    integer(i4b)                     :: iFile=1                    ! index of current forcing file from forcing file list
    integer(i4b)                     :: forcingStep=integerMissing ! index of current time step in current forcing file
    real(dp),allocatable             :: zSoilReverseSign(:)        ! height at bottom of each soil layer, negative downwards (m)
    real(dp),dimension(12)           :: greenVegFrac_monthly       ! fraction of green vegetation in each month (0-1)
    logical(lgt),parameter           :: overwriteRSMIN=.false.     ! flag to overwrite RSMIN
    real(dp)                         :: notUsed_canopyDepth        ! NOT USED: canopy depth (m)
    real(dp)                         :: notUsed_exposedVAI         ! NOT USED: exposed vegetation area index (m2 m-2)
    ! error control
    integer(i4b)                     :: err=0                      ! error code
    character(len=1024)              :: message=''                 ! error message
    ! output control 
    integer(i4b)                     :: iFreq                      ! index for looping through output files
    logical(lgt)                     :: statForc_mask(maxvarForc)  ! mask defining forc stats
    logical(lgt)                     :: statProg_mask(maxvarProg)  ! mask defining prog stats
    logical(lgt)                     :: statDiag_mask(maxvarDiag)  ! mask defining diag stats
    logical(lgt)                     :: statFlux_mask(maxvarFlux)  ! mask defining flux stats
    logical(lgt)                     :: statIndx_mask(maxvarIndx)  ! mask defining indx stats
    logical(lgt)                     :: statBvar_mask(maxvarBvar)  ! mask defining bvar stats
    integer(i4b),allocatable         :: forcChild_map(:)           ! index of the child data structure: stats forc
    integer(i4b),allocatable         :: progChild_map(:)           ! index of the child data structure: stats prog
    integer(i4b),allocatable         :: diagChild_map(:)           ! index of the child data structure: stats diag
    integer(i4b),allocatable         :: fluxChild_map(:)           ! index of the child data structure: stats flux
    integer(i4b),allocatable         :: indxChild_map(:)           ! index of the child data structure: stats indx
    integer(i4b),allocatable         :: bvarChild_map(:)           ! index of the child data structure: stats bvar
    type(extended_info),allocatable  :: statForc_meta(:)           ! child metadata for stats 
    type(extended_info),allocatable  :: statProg_meta(:)           ! child metadata for stats 
    type(extended_info),allocatable  :: statDiag_meta(:)           ! child metadata for stats 
    type(extended_info),allocatable  :: statFlux_meta(:)           ! child metadata for stats 
    type(extended_info),allocatable  :: statIndx_meta(:)           ! child metadata for stats 
    type(extended_info),allocatable  :: statBvar_meta(:)           ! child metadata for stats 
    ! stuff for restart file
    character(len=256)               :: timeString                 ! protion of restart file name that contains the write-out time
    character(len=256)               :: restartFile                ! restart file name
    character(len=256)               :: attrFile                   ! attributes file name
    ! parallelize the model run
    integer(i4b)                     :: startGRU                   ! index of the starting GRU for parallelization run
    integer(i4b)                     :: checkHRU                   ! index of the HRU for a single HRU run
    integer(i4b)                     :: fileGRU                    ! number of GRUs in the input file
    integer(i4b)                     :: fileHRU                    ! number of HRUs in the input file
    integer(i4b)                     :: iRunMode                   ! define the current running mode
    character(len=128)               :: fmtGruOutput               ! a format string used to write start and end GRU in output file names

    contains

end module bmi_summa
