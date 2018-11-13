SELECT

  main.object_id, main.ra, main.dec, main.parent_id, main.deblend_nchild,

  shape_detradius(array[main.gshape_sdss_psf_11, main.gshape_sdss_psf_22, main.gshape_sdss_psf_12]) as gPSF_pix,
  shape_detradius(array[main.rshape_sdss_psf_11, main.rshape_sdss_psf_22, main.rshape_sdss_psf_12]) as rPSF_pix,
  shape_detradius(array[main.ishape_sdss_psf_11, main.ishape_sdss_psf_22, main.ishape_sdss_psf_12]) as iPSF_pix,
  shape_detradius(array[main.zshape_sdss_psf_11, main.zshape_sdss_psf_22, main.zshape_sdss_psf_12]) as zPSF_pix,
  shape_detradius(array[main.yshape_sdss_psf_11, main.yshape_sdss_psf_22, main.yshape_sdss_psf_12]) as yPSF_pix,

  shape_detradius(array[main.gshape_sdss_11, main.gshape_sdss_22, main.gshape_sdss_12]) as gSize_pix,
  shape_detradius(array[main.rshape_sdss_11, main.rshape_sdss_22, main.rshape_sdss_12]) as rSize_pix,
  shape_detradius(array[main.ishape_sdss_11, main.ishape_sdss_22, main.ishape_sdss_12]) as iSize_pix,
  shape_detradius(array[main.zshape_sdss_11, main.zshape_sdss_22, main.zshape_sdss_12]) as zSize_pix,
  shape_detradius(array[main.yshape_sdss_11, main.yshape_sdss_22, main.yshape_sdss_12]) as ySize_pix,

  main.a_g, main.a_r, main.a_i, main.a_z, main.a_y,

  main.gmag_psf, main.rmag_psf, main.imag_psf, main.zmag_psf, main.ymag_psf,
  main.gmag_psf_err, main.rmag_psf_err, main.imag_psf_err, main.zmag_psf_err, main.ymag_psf_err,
  main.gflux_psf_flags, main.rflux_psf_flags, main.iflux_psf_flags, main.zflux_psf_flags, main.yflux_psf_flags,

  main.gmag_kron, main.rmag_kron, main.imag_kron, main.zmag_kron, main.ymag_kron,
  main.gmag_kron_err, main.rmag_kron_err, main.imag_kron_err, main.zmag_kron_err, main.ymag_kron_err,
  main.gflux_kron_flags, main.rflux_kron_flags, main.iflux_kron_flags, main.zflux_kron_flags, main.yflux_kron_flags,

  main.gcmodel_mag, main.rcmodel_mag, main.icmodel_mag, main.zcmodel_mag, main.ycmodel_mag,
  main.gcmodel_mag_err, main.rcmodel_mag_err, main.icmodel_mag_err, main.zcmodel_mag_err, main.ycmodel_mag_err,
  main.gcmodel_flux_flags, main.rcmodel_flux_flags, main.icmodel_flux_flags, main.zcmodel_flux_flags, main.ycmodel_flux_flags,

  main.gcountinputs, main.rcountinputs, main.icountinputs, main.zcountinputs, main.ycountinputs,

  main.gvariance, main.rvariance, main.ivariance, main.zvariance, main.yvariance,

  main.gclassification_extendedness, main.rclassification_extendedness, main.iclassification_extendedness, main.zclassification_extendedness, main.yclassification_extendedness,

  main.gflags_pixel_edge, main.rflags_pixel_edge, main.iflags_pixel_edge, main.zflags_pixel_edge, main.yflags_pixel_edge,
  main.gflags_pixel_interpolated_any, main.rflags_pixel_interpolated_any, main.iflags_pixel_interpolated_any, main.zflags_pixel_interpolated_any, main.yflags_pixel_interpolated_any,
  main.gflags_pixel_interpolated_center, main.rflags_pixel_interpolated_center, main.iflags_pixel_interpolated_center, main.zflags_pixel_interpolated_center, main.yflags_pixel_interpolated_center,
  main.gflags_pixel_saturated_any, main.rflags_pixel_saturated_any, main.iflags_pixel_saturated_any, main.zflags_pixel_saturated_any, main.yflags_pixel_saturated_any,
  main.gflags_pixel_saturated_center, main.rflags_pixel_saturated_center, main.iflags_pixel_saturated_center, main.zflags_pixel_saturated_center, main.yflags_pixel_saturated_center,
  main.gflags_pixel_cr_any, main.rflags_pixel_cr_any, main.iflags_pixel_cr_any, main.zflags_pixel_cr_any, main.yflags_pixel_cr_any,
  main.gflags_pixel_cr_center, main.rflags_pixel_cr_center, main.iflags_pixel_cr_center, main.zflags_pixel_cr_center, main.yflags_pixel_cr_center,
  main.gflags_pixel_bad, main.rflags_pixel_bad, main.iflags_pixel_bad, main.zflags_pixel_bad, main.yflags_pixel_bad,
  main.gflags_pixel_suspect_any, main.rflags_pixel_suspect_any, main.iflags_pixel_suspect_any, main.zflags_pixel_suspect_any, main.yflags_pixel_suspect_any,
  main.gflags_pixel_suspect_center, main.rflags_pixel_suspect_center, main.iflags_pixel_suspect_center, main.zflags_pixel_suspect_center, main.yflags_pixel_suspect_center,
  main.gflags_pixel_offimage, main.rflags_pixel_offimage, main.iflags_pixel_offimage, main.zflags_pixel_offimage, main.yflags_pixel_offimage,
  main.gflags_pixel_bright_object_center, main.rflags_pixel_bright_object_center, main.iflags_pixel_bright_object_center, main.zflags_pixel_bright_object_center, main.yflags_pixel_bright_object_center,
  main.gflags_pixel_clipped_any, main.rflags_pixel_clipped_any, main.iflags_pixel_clipped_any, main.zflags_pixel_clipped_any, main.yflags_pixel_clipped_any,
  main.gflags_pixel_bright_object_any, main.rflags_pixel_bright_object_any, main.iflags_pixel_bright_object_any, main.zflags_pixel_bright_object_any, main.yflags_pixel_bright_object_any

FROM
  s16a_wide.forced AS main

WHERE
  (main.detect_is_tract_inner = 't'
  AND main.detect_is_patch_inner = 't')
  AND (main.merge_peak_sky = 'f')
  AND (main.gmag_psf < 20 OR main.rmag_psf < 20 OR main.imag_psf < 20 OR main.zmag_psf < 20 OR main.ymag_psf < 20 )
