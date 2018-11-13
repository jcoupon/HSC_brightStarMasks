WITH
    {user_catalog}
    ,
    match AS (
        SELECT
            object_id,
            earth_distance(coord, ll_to_earth({dec}, {ra})) AS match_distance,
            user_catalog.*
        FROM
            user_catalog JOIN {rerun}.forced
                ON coneSearch(coord, {ra}, {dec}, {accuracy})
    )
SELECT
   match.*,
   object_id, ra, dec, parent_id, deblend_nchild,
   shape_detradius(array[gshape_sdss_psf_11, gshape_sdss_psf_22, gshape_sdss_psf_12])*2.0*sqrt(2.0) as gPSF_FWHM,
   shape_detradius(array[rshape_sdss_psf_11, rshape_sdss_psf_22, rshape_sdss_psf_12])*2.0*sqrt(2.0) as rPSF_FWHM,
   shape_detradius(array[ishape_sdss_psf_11, ishape_sdss_psf_22, ishape_sdss_psf_12])*2.0*sqrt(2.0) as iPSF_FWHM,
   shape_detradius(array[zshape_sdss_psf_11, zshape_sdss_psf_22, zshape_sdss_psf_12])*2.0*sqrt(2.0) as zPSF_FWHM,
   shape_detradius(array[yshape_sdss_psf_11, yshape_sdss_psf_22, yshape_sdss_psf_12])*2.0*sqrt(2.0) as yPSF_FWHM,

   gmag_psf, rmag_psf, imag_psf, zmag_psf, ymag_psf,
   gmag_psf_err, rmag_psf_err, imag_psf_err, zmag_psf_err, ymag_psf_err,
   gflux_psf_flags, rflux_psf_flags, iflux_psf_flags, zflux_psf_flags, yflux_psf_flags,

   gcmodel_mag, rcmodel_mag, icmodel_mag, zcmodel_mag, ycmodel_mag,
   gcmodel_mag_err, rcmodel_mag_err, icmodel_mag_err, zcmodel_mag_err, ycmodel_mag_err,
   gcmodel_flux_flags, rcmodel_flux_flags, icmodel_flux_flags, zcmodel_flux_flags, ycmodel_flux_flags,

   gflags_pixel_edge, rflags_pixel_edge, iflags_pixel_edge, zflags_pixel_edge, yflags_pixel_edge,
   gflags_pixel_interpolated_any, rflags_pixel_interpolated_any, iflags_pixel_interpolated_any, zflags_pixel_interpolated_any, yflags_pixel_interpolated_any,
   gflags_pixel_interpolated_center, rflags_pixel_interpolated_center, iflags_pixel_interpolated_center, zflags_pixel_interpolated_center, yflags_pixel_interpolated_center,
   gflags_pixel_saturated_any, rflags_pixel_saturated_any, iflags_pixel_saturated_any, zflags_pixel_saturated_any, yflags_pixel_saturated_any,
   gflags_pixel_saturated_center, rflags_pixel_saturated_center, iflags_pixel_saturated_center, zflags_pixel_saturated_center, yflags_pixel_saturated_center,
   gflags_pixel_cr_any, rflags_pixel_cr_any, iflags_pixel_cr_any, zflags_pixel_cr_any, yflags_pixel_cr_any,
   gflags_pixel_cr_center, rflags_pixel_cr_center, iflags_pixel_cr_center, zflags_pixel_cr_center, yflags_pixel_cr_center,
   gflags_pixel_bad, rflags_pixel_bad, iflags_pixel_bad, zflags_pixel_bad, yflags_pixel_bad,
   gflags_pixel_suspect_any, rflags_pixel_suspect_any, iflags_pixel_suspect_any, zflags_pixel_suspect_any, yflags_pixel_suspect_any,
   gflags_pixel_suspect_center, rflags_pixel_suspect_center, iflags_pixel_suspect_center, zflags_pixel_suspect_center, yflags_pixel_suspect_center,
   gflags_pixel_offimage, rflags_pixel_offimage, iflags_pixel_offimage, zflags_pixel_offimage, yflags_pixel_offimage,
   gflags_pixel_bright_object_center, rflags_pixel_bright_object_center, iflags_pixel_bright_object_center, zflags_pixel_bright_object_center, yflags_pixel_bright_object_center,
   gflags_pixel_clipped_any, rflags_pixel_clipped_any, iflags_pixel_clipped_any, zflags_pixel_clipped_any, yflags_pixel_clipped_any,
   gflags_pixel_bright_object_any, rflags_pixel_bright_object_any, iflags_pixel_bright_object_any, zflags_pixel_bright_object_any, yflags_pixel_bright_object_any

   {columns}
FROM
    match LEFT JOIN {rerun}.forced USING(object_id)
WHERE
   (detect_is_tract_inner = 't' AND detect_is_patch_inner = 't' AND merge_peak_sky = 'f')
