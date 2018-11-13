WITH
    {user_catalog}
    ,
    match AS (
        SELECT
            object_id,
            earth_distance(ll_to_earth(dec, ra), ll_to_earth({dec}, {ra})) AS match_distance,
            user_catalog.*
        FROM
            user_catalog JOIN {rerun}.random
                ON coneSearch(ll_to_earth(dec, ra), {ra}, {dec}, {accuracy})
    )
SELECT
   match.*,
   object_id, ra, dec

   {columns}
FROM
    match LEFT JOIN {rerun}.random USING(object_id)
WHERE
   (idetect_is_tract_inner = 't' AND idetect_is_patch_inner = 't')
