MODULE SUMMALIB

    USE nrtype                                                  ! variable types, etc.
    USE netcdf                                                  ! netcdf libraries
    USE,intrinsic :: ieee_arithmetic                            ! IEEE arithmetic (obviously)
    ! provide access to subroutines and functions
    USE summaFileManager,only:summa_SetDirsUndPhiles            ! sets directories and filenames
    USE module_sf_noahmplsm,only:read_mp_veg_parameters         ! module to read NOAH vegetation tables
    USE module_sf_noahmplsm,only:isWater                        ! parameter for water land cover type
    USE nr_utility_module,only:arth                             ! get a sequence of numbers
    USE nr_utility_module,only:indexx                           ! sort vectors in ascending order
    USE ascii_util_module,only:file_open                        ! open ascii file
    USE ascii_util_module,only:get_vlines                       ! read a vector of non-comment lines from an ASCII file
    USE ascii_util_module,only:split_line                       ! extract the list of variable names from the character string
    use time_utils_module,only:elapsedSec                       ! calculate the elapsed time
    USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures
    USE allocspace_module,only:allocLocal                       ! module to allocate space for local data structures
    USE childStruc_module,only:childStruc                       ! module to create a child data structure
    USE mDecisions_module,only:mDecisions                       ! module to read model decisions
    USE popMetadat_module,only:popMetadat                       ! module to populate metadata structures
    USE flxMapping_module,only:flxMapping                       ! module to map fluxes to states
    USE checkStruc_module,only:checkStruc                       ! module to check metadata structures
    USE def_output_module,only:def_output                       ! module to define model output
    USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
    USE read_attrb_module,only:read_dimension                   ! module to read dimensions of GRU and HRU
    USE read_attrb_module,only:read_attrb                       ! module to read local attributes
    USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
    USE paramCheck_module,only:paramCheck                       ! module to check consistency of model parameters
    USE check_icond_module,only:check_icond                     ! module to check initial conditions
    USE read_icond_module,only:read_icond                       ! module to read initial conditions
    USE read_icond_module,only:read_icond_nlayers               ! module to read initial conditions
    USE pOverwrite_module,only:pOverwrite                       ! module to overwrite default parameter values with info from the Noah tables
    USE read_param_module,only:read_param                       ! module to read model parameter sets
    USE ConvE2Temp_module,only:E2T_lookup                       ! module to calculate a look-up table for the temperature-enthalpy conversion
    USE var_derive_module,only:calcHeight                       ! module to calculate height at layer interfaces and layer mid-point
    USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
    USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
    USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
    USE var_derive_module,only:fracFuture                       ! module to calculate the fraction of runoff in future time steps (time delay histogram)
    USE read_force_module,only:read_force                       ! module to read model forcing data
    USE modelwrite_module,only:writeParm,writeTime              ! module to write model attributes and parameters
    USE modelwrite_module,only:writeData,writeBasin             ! module to write model output
    USE modelwrite_module,only:writeRestart                     ! module to write model Restart
    USE vegPhenlgy_module,only:vegPhenlgy                       ! module to compute vegetation phenology
    USE run_oneGRU_module,only:run_oneGRU                       ! module to run for one GRU
    USE groundwatr_module,only:groundwatr                       ! module to simulate regional groundwater balance
    USE qTimeDelay_module,only:qOverland                        ! module to route water through an "unresolved" river network
    USE netcdf_util_module,only:nc_file_close                   ! module to handle netcdf stuff for inputs and outputs
    ! provide access to file paths
    USE summaFileManager,only:SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
    USE summaFileManager,only:MODEL_INITCOND                    ! name of model initial conditions file
    USE summaFileManager,only:LOCAL_ATTRIBUTES                  ! name of model initial attributes file
    USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
    USE summaFileManager,only:LOCALPARAM_INFO,BASINPARAM_INFO   ! files defining the default values and constraints for model parameters
    ! provide access to the derived types to define the data structures
    USE data_types,only:&
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
    USE data_types,only:extended_info          ! extended metadata structure
    ! provide access to runtime options
    USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU
    ! provide access to metadata structures
    USE globalData,only:time_meta,forc_meta,attr_meta,type_meta ! metadata structures
    USE globalData,only:prog_meta,diag_meta,flux_meta           ! metadata structures
    USE globalData,only:mpar_meta,indx_meta                     ! metadata structures
    USE globalData,only:bpar_meta,bvar_meta                     ! metadata structures
    USE globalData,only:averageFlux_meta                        ! metadata for time-step average fluxes
    USE globalData,only:model_decisions                         ! model decision structure
    ! provide access to global data
    USE globalData,only:dNaN                                    ! double precision NaN
    USE globalData,only:refTime                                 ! reference time
    USE globalData,only:startTime                               ! start time
    USE globalData,only:finshTime                               ! end time
    USE globalData,only:doJacobian                              ! flag to compute the Jacobian
    USE globalData,only:gru_struc                               ! gru-hru mapping structures
    USE globalData,only:localParFallback                        ! local column default parameters
    USE globalData,only:basinParFallback                        ! basin-average default parameters
    USE globalData,only:structInfo                              ! information on the data structures
    USE globalData,only:numtim                                  ! number of time steps
    USE globalData,only:urbanVegCategory                        ! vegetation category for urban areas
    USE globalData,only:greenVegFrac_monthly                    ! fraction of green vegetation in each month (0-1)
    USE globalData,only:globalPrintFlag                         ! global print flag
    USE globalData,only:integerMissing                          ! missing integer value
    USE globalData,only:realMissing                             ! missing double precision value
    USE globalData,only:yes,no                                  ! .true. and .false.
    ! provide access to Noah-MP parameters
    USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM                    ! 2-d tables for stem area index and leaf area index (vegType,month)
    USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB                      ! height at the top and bottom of vegetation (vegType)
    USE var_lookup,only:maxvarTime                              ! size of variable vectors
    USE var_lookup,only:maxvarForc,maxvarProg,maxvarDiag        ! size of variable vectors
    USE var_lookup,only:maxvarFlux,maxvarIndx,maxvarBvar        ! size of variable vectors
    ! provide access to the named variables that describe elements of parent model structures
    USE var_lookup,only:iLookTIME,iLookFORCE                    ! look-up values for time and forcing data structures
    USE var_lookup,only:iLookTYPE                               ! look-up values for classification of veg, soils etc.
    USE var_lookup,only:iLookATTR                               ! look-up values for local attributes
    USE var_lookup,only:iLookPARAM                              ! look-up values for local column model parameters
    USE var_lookup,only:iLookINDEX                              ! look-up values for local column index variables
    USE var_lookup,only:iLookPROG                               ! look-up values for local column model prognostic (state) variables
    USE var_lookup,only:iLookDIAG                               ! look-up values for local column model diagnostic variables
    USE var_lookup,only:iLookFLUX                               ! look-up values for local column model fluxes
    USE var_lookup,only:iLookBVAR                               ! look-up values for basin-average model variables
    USE var_lookup,only:iLookBPAR                               ! look-up values for basin-average model parameters
    USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions
    USE var_lookup,only:iLookVarType                            ! look-up values for variable type structure
    USE var_lookup,only:iLookFreq                               ! look-up values for model output frequency
    ! provide access to the named variables that describe elements of child  model structures
    USE var_lookup,only:childFLUX_MEAN                          ! look-up values for timestep-average model fluxes
    ! provide access to the named variables that describe model decisions
    USE mDecisions_module,only:  &                              ! look-up values for method used to compute derivative
     numerical,   & ! numerical solution
     analytical     ! analytical solution
    USE mDecisions_module,only:&                                ! look-up values for LAI decisions
     monthlyTable,& ! LAI/SAI taken directly from a monthly table for different vegetation classes
     specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
    USE mDecisions_module,only:&                                ! look-up values for the choice of method for the spatial representation of groundwater
     localColumn, & ! separate groundwater representation in each local soil column
     singleBasin    ! single groundwater store over the entire basin
    USE mDecisions_module,only:&
      sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
      rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
    USE output_stats,only:calcStats                             ! module for compiling output statistics
    USE var_lookup,only:maxvarFreq                              ! maximum # of output files
    USE globalData,only:ncid                                    ! file id of netcdf output file

    implicit none

CONTAINS

    FUNCTION get_summa_time() RESULT(ret) bind(c, name="get_summa_time")
    END FUNCTION get_summa_time


    FUNCTION get_end_time() RESULT(ret) bind(c, name="get_summa_end_time")
    END FUNCTION get_end_time


    FUNCTION get_time_step() RESULT(ret) bind(c, name="get_summa_time_step")
    END FUNCTION get_time_step


    FUNCTION convert_summa_time(idate) RESULT(ret)
    END FUNCTION convert_summa_time


    FUNCTION get_num_output_fields() RESULT(ret) bind(c, name="get_num_ovars")
    END FUNCTION get_num_output_fields


    SUBROUTINE get_output_name(index, dest) bind(c, name="get_ovar_name")
    END SUBROUTINE get_output_name


    SUBROUTINE get_output_units(index, dest) bind(c, name="get_ovar_units")
    END SUBROUTINE get_output_units


    FUNCTION get_num_subbasins() RESULT(ret) bind(c, name="get_num_basins")
    END FUNCTION get_num_subbasins


    SUBROUTINE get_latlons(targetlatarr, targetlonarr) bind(c, name="get_latlons")
    END SUBROUTINE get_latlons


    SUBROUTINE get_basin_field(index, targetarr) bind(c, name="get_ovar_values")
    END SUBROUTINE get_basin_field


    FUNCTION initialize(dir, iseq) RESULT(istat) bind(c, name="init_summa")
    END FUNCTION initialize


    FUNCTION update() RESULT(istat) bind(c, name="update_summa")
    END FUNCTION update


    FUNCTION finalize() RESULT(istat) bind(c, name="finalize_summa")
    END FUNCTION finalize


END MODULE SUMMALIB
