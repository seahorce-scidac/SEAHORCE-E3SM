# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Modules"

list ( APPEND _files
       Modules/mod_arrays.F
       Modules/mod_average.F
       Modules/mod_bbl.F
       Modules/mod_behavior.F
       Modules/mod_biology.F
       Modules/mod_boundary.F
       Modules/mod_clima.F
       Modules/mod_coupler.F
       Modules/mod_coupling.F
       Modules/mod_diags.F
       Modules/mod_eclight.F
       Modules/mod_eoscoef.F
       Modules/mod_floats.F
       Modules/mod_forces.F
       Modules/mod_fourdvar.F
       Modules/mod_grid.F
       Modules/mod_iounits.F
       Modules/mod_kinds.F
       Modules/mod_mixing.F
       Modules/mod_ncparam.F
       Modules/mod_nesting.F
       Modules/mod_netcdf.F
       Modules/mod_ocean.F
       Modules/mod_parallel.F
       Modules/mod_param.F
       Modules/mod_scalars.F
       Modules/mod_sedbed.F
       Modules/mod_sediment.F
       Modules/mod_sources.F
       Modules/mod_stepping.F
       Modules/mod_storage.F
       Modules/mod_strings.F
       Modules/mod_tides.F
)

set ( ROMS_Modules_files
      ${_files}
      PARENT_SCOPE
)
