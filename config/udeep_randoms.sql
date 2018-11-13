

SELECT main.ra, main.dec,
  shape_detradius(array[main.gshape_sdss_psf_xx, main.gshape_sdss_psf_yy, main.gshape_sdss_psf_xy]) as gPSF_pix,
  shape_detradius(array[main.rshape_sdss_psf_xx, main.rshape_sdss_psf_yy, main.rshape_sdss_psf_xy]) as rPSF_pix,
  shape_detradius(array[main.ishape_sdss_psf_xx, main.ishape_sdss_psf_yy, main.ishape_sdss_psf_xy]) as iPSF_pix,
  shape_detradius(array[main.zshape_sdss_psf_xx, main.zshape_sdss_psf_yy, main.zshape_sdss_psf_xy]) as zPSF_pix,
  shape_detradius(array[main.yshape_sdss_psf_xx, main.yshape_sdss_psf_yy, main.yshape_sdss_psf_xy]) as yPSF_pix,

  main.gcountinputs, main.rcountinputs, main.icountinputs, main.zcountinputs, main.ycountinputs,

  main.gpix_variance, main.rpix_variance, main.ipix_variance, main.zpix_variance, main.ypix_variance,

  main.gflags_pixel_edge, main.rflags_pixel_edge, main.iflags_pixel_edge, main.zflags_pixel_edge, main.yflags_pixel_edge,
  main.gflags_pixel_interpolated_center, main.rflags_pixel_interpolated_center, main.iflags_pixel_interpolated_center, main.zflags_pixel_interpolated_center, main.yflags_pixel_interpolated_center,
  main.gflags_pixel_saturated_center, main.rflags_pixel_saturated_center, main.iflags_pixel_saturated_center, main.zflags_pixel_saturated_center, main.yflags_pixel_saturated_center,
  main.gflags_pixel_cr_center, main.rflags_pixel_cr_center, main.iflags_pixel_cr_center, main.zflags_pixel_cr_center, main.yflags_pixel_cr_center,
  main.gflags_pixel_bad, main.rflags_pixel_bad, main.iflags_pixel_bad, main.zflags_pixel_bad, main.yflags_pixel_bad,
  main.gflags_pixel_suspect_center, main.rflags_pixel_suspect_center, main.iflags_pixel_suspect_center, main.zflags_pixel_suspect_center, main.yflags_pixel_suspect_center,
  main.gflags_pixel_offimage, main.rflags_pixel_offimage, main.iflags_pixel_offimage, main.zflags_pixel_offimage, main.yflags_pixel_offimage,
  main.gflags_pixel_bright_object_center, main.rflags_pixel_bright_object_center, main.iflags_pixel_bright_object_center, main.zflags_pixel_bright_object_center, main.yflags_pixel_bright_object_center

from s16a_udeep_random.random AS main

where (main.idetect_is_tract_inner = 't' AND main.idetect_is_patch_inner = 't') AND main.iadjust_density < 0.01
