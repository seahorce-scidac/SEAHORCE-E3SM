# svn $Id: CMakeLists.txt 1054 2021-03-06 19:47:12Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::: David Robertson :::
# Copyright (c) 2002-2021 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Source code list for sub-directory "ROMS/Nonlinear"

list( APPEND _files
      Nonlinear/bbl.F
      Nonlinear/bc_2d.F
      Nonlinear/bc_3d.F
      Nonlinear/bc_4d.F
      Nonlinear/bc_bry2d.F
      Nonlinear/bc_bry3d.F
      Nonlinear/bulk_flux.F
      Nonlinear/bvf_mix.F
      Nonlinear/conv_2d.F
      Nonlinear/conv_3d.F
      Nonlinear/conv_bry2d.F
      Nonlinear/conv_bry3d.F
      Nonlinear/diag.F
      Nonlinear/exchange_2d.F
      Nonlinear/exchange_3d.F
      Nonlinear/exchange_4d.F
      Nonlinear/forcing.F
      Nonlinear/frc_adjust.F
      Nonlinear/get_data.F
      Nonlinear/get_idata.F
      Nonlinear/gls_corstep.F
      Nonlinear/gls_prestep.F
      Nonlinear/hmixing.F
      Nonlinear/ini_fields.F
      Nonlinear/initial.F
      Nonlinear/interp_floats.F
      Nonlinear/lmd_bkpp.F
      Nonlinear/lmd_skpp.F
      Nonlinear/lmd_swfrac.F
      Nonlinear/lmd_vmix.F
      Nonlinear/main2d.F
      Nonlinear/main3d.F
      Nonlinear/mpdata_adiff.F
      Nonlinear/my25_corstep.F
      Nonlinear/my25_prestep.F
      Nonlinear/nesting.F
      Nonlinear/obc_adjust.F
      Nonlinear/obc_volcons.F
      Nonlinear/omega.F
      Nonlinear/output.F
      Nonlinear/pre_step3d.F
      Nonlinear/prsgrd.F
      # Nonlinear/radiation_stress.F
      Nonlinear/rho_eos.F
      Nonlinear/rhs3d.F
      Nonlinear/set_avg.F
      Nonlinear/set_data.F
      Nonlinear/set_depth.F
      Nonlinear/set_massflux.F
      Nonlinear/set_tides.F
      Nonlinear/set_vbc.F
      Nonlinear/set_zeta.F
      Nonlinear/step2d.F
      Nonlinear/step3d_t.F
      Nonlinear/step3d_uv.F
      Nonlinear/step_floats.F
      Nonlinear/t3dbc_im.F
      Nonlinear/t3dmix.F
      Nonlinear/tkebc_im.F
      Nonlinear/u2dbc_im.F
      Nonlinear/u3dbc_im.F
      Nonlinear/uv3dmix.F
      Nonlinear/v2dbc_im.F
      Nonlinear/v3dbc_im.F
      Nonlinear/vwalk_floats.F
      Nonlinear/wetdry.F
      Nonlinear/wvelocity.F
      Nonlinear/zetabc.F
)

set ( ROMS_Nonlinear_files
      ${_files}
      PARENT_SCOPE
)
