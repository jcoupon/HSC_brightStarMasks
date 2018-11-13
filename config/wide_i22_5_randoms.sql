SELECT
  main.object_id, main.ra, main.dec,
  main.icountinputs,
  main.ipix_variance,
  main.iflags_pixel_bright_object_center,
  main.iflags_pixel_bright_object_any

FROM
  s16a_wide_random.random AS main

WHERE
   NOT main.iflags_pixel_edge
   AND main.idetect_is_tract_inner = 't'
   AND main.idetect_is_patch_inner = 't'
   AND iadjust_density < 0.1
