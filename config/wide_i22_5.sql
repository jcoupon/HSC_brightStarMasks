SELECT
  main.object_id, main.ra, main.dec,
  main.imag_kron,
  main.imag_kron_err,
  main.icmodel_mag,
  main.icmodel_mag_err,
  main.icmodel_flux_flags,
  main.icountinputs,
  main.ivariance,
  main.iclassification_extendedness,
  main.iflags_pixel_bright_object_center,
  main.iflags_pixel_bright_object_any,
  main.a_i
FROM
  s16a_wide.forced AS main

WHERE
   main.icmodel_mag-main.a_i < 22.5
   AND NOT main.icentroid_sdss_flags
   AND NOT main.iflags_pixel_edge
   AND NOT main.iflags_pixel_interpolated_any
   AND NOT main.iflags_pixel_saturated_any
   AND NOT main.iflags_pixel_cr_any
   AND NOT main.iflags_pixel_bad
   AND NOT main.merge_peak_sky
   AND main.detect_is_primary
