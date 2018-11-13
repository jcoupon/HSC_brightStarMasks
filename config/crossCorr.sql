WITH
    user_catalog("user_ra","user_dec","user_G_Gaia","user_r_Czakon_deg","user_ipsf_FWHM") AS (VALUES('3.3027089999999998e+02'::double precision,'6.0471999999999992e-01'::double precision,'4.9366323736854749e+00'::double precision,'1.9412253526607853e-01'::double precision,'6.7521903747023015e-01'::double precision),('3.3631925999999993e+02','1.3773999999999997e+00','4.8120494611860787e+00','2.0784184553592594e-01','6.4896212254553021e-01'),('3.4836619999999996e+01','-2.9776499999999997e+00','4.0034374520556852e+00','3.2505925992598350e-01','8.4733090373290065e-01'),('2.2275446999999997e+02','-2.2991499999999996e+00','4.6226870502929733e+00','2.3065320582364168e-01','6.7790665915743364e-01'),('2.1488533999999999e+02','-2.2655199999999995e+00','4.8742613656997431e+00','2.0086921023220544e-01','8.2007598177367302e-01'),('2.1705056999999996e+02','-2.2279499999999999e+00','4.6770277896545362e+00','2.2385318495482701e-01','7.9423837526305752e-01'),('2.2129877999999997e+02','-1.4175299999999997e+00','4.8933327126910164e+00','1.9878075519890001e-01','5.8270820841459881e-01'),('1.2968931999999998e+02','3.3414399999999995e+00','4.1127097032951063e+00','3.0588499434033062e-01','0.0000000000000000e+00'),('1.3080614999999997e+02','3.3986599999999996e+00','4.2920494802595659e+00','2.7689906976071027e-01','0.0000000000000000e+00'),('1.3649318999999997e+02','5.0923199999999991e+00','4.6627097148151444e+00','2.2562449974743284e-01','5.2242342382389539e-01'),('2.4852574999999996e+02','4.2437029999999993e+01','4.1878641606108831e+00','2.9337561500812098e-01','5.3184059614623280e-01'),('2.3816892999999996e+02','4.2451519999999995e+01','4.3842744233879820e+00','2.6311010825391490e-01','6.5529263972430329e-01'),('2.3865770999999998e+02','4.3138569999999994e+01','4.2233326363970711e+00','2.8765594920305670e-01','6.8286487093129700e-01'),('2.1337085999999996e+02','5.1789969999999997e+01','4.5083940574355399e+00','2.4566235545671092e-01','0.0000000000000000e+00'))
    ,
    match AS (
        SELECT
            object_id,
            earth_distance(coord, ll_to_earth("user_dec", "user_ra")) AS match_distance,
            user_catalog.*
        FROM
            user_catalog JOIN "s16a_wide".forced
                ON coneSearch(coord, "user_ra", "user_dec", 1200.0)
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

   
FROM
    match LEFT JOIN "s16a_wide".forced USING(object_id)
WHERE
   (detect_is_tract_inner = 't' AND detect_is_patch_inner = 't' AND merge_peak_sky = 'f')
