#! /bin/bash
set -o nounset


source scripts/config.sh

CASJOB='java -jar /Users/coupon/local/bin/casjobs.jar'


# HSC field limits (ra_min, ra_max, dec_min, dec_max)

# coordinates
# fall1 21h56m 24h00m -02d00m 08d00m
# fall2 00h00m 2h44m -02d00m 08d00m
# fall3 01h46 02h44m -08d00m -02d00m
# Spring 08h26m 15h04 -03d00m 06d00
# Northern 13h16m 16h44m 41.5d00m 45d00m
# Elais center: 16h10 55d00 r = 2.5deg
# AEGIS center: 14h18m 52d30m r = 2.0deg
# pipeTest 21h08m  21h36m  -03d00m 04d00m

# decimal degrees
# fall1 329.00 360.0 -02.00 08.00
# fall2 0.0 41.0 -02.00 08.00
# fall3 26.50 41.00 -08.00 -02.00
# Spring 126.5 226.0 -03.00 06.00
# Northern 199.0 251.00 41.5 45.00
# Elais center: 242.5 55.0 r = 2.5deg
# AEGIS center: 214.50 52.50 r = 2.0deg
# pipeTest 317 324 -3.0 4.0



main() {

   # getGaia
   # getIGSL
   # getPickles
   # getHSC
   # getCFHTLenS
   # getPanSTARRS
   getSDSS

   return
}

# ----------------------------------------------------------------------------------- #
# functions
# ----------------------------------------------------------------------------------- #


function getGaia
{

   # http://www.robertmartinayers.org/tools/coordinates.html
   # http://gea.esac.esa.int/archive/

   # Fall
   # SELECT * FROM gaiadr1.gaia_source
   # WHERE
   # CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',344.5,3,31,10))=1
   # OR
   # CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',20.5,3,41,10))=1
   # OR
   # CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',33.75,-5,14.5,6))=1

   # Spring
   # SELECT * FROM gaiadr1.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',176.25,1.5,99.5,9))=1

   # Northern
   # SELECT * FROM gaiadr1.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',225,43.25,52,3.5))=1

   # Elais
   # SELECT * FROM gaiadr1.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),CIRCLE('ICRS',242.5,54.0,2.5))=1

   # AEGIS
   # SELECT * FROM gaiadr1.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),CIRCLE('ICRS',214.50,52.50,2.0))=1

   $STILTS  tcatn nin=5 in1=$DATADIR/Gaia/Fall.fits in2=$DATADIR/Gaia/Spring.fits in3=$DATADIR/Gaia/Northern.fits in4=$DATADIR/Gaia/Elais.fits in5=$DATADIR/Gaia/Aegis.fits out=$DATADIR/Gaia/Gaia_HSC_footprint.fits

   # pipeTest
   # SELECT * FROM gaiadr1.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),BOX('ICRS',320.5,0.5,7.0,7.0))=1


   return
}


function getIGSL
{

   # http://www.robertmartinayers.org/tools/coordinates.html
   # http://gea.esac.esa.int/archive/

   # Fall
   # SELECT * FROM public.igsl_source WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',344.5,3,31,10))=1
   # SELECT * FROM public.igsl_source WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',13.75,3,27.5,10))=1
   # SELECT * FROM public.igsl_source WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',34.25,3,13.5,10))=1
   # SELECT * FROM public.igsl_source WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',33.75,-5,14.5,6))=1

   # Spring
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',138.25,1.5,23.5,9))=1
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',161.25,1.5,22.5,9))=1
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',183.75,1.5,22.5,9))=1
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',206.25,1.5,22.5,9))=1
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',221.75,1.5,8.5,9))=1

   # Northern
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),BOX('ICRS',225,43.25,52,3.5))=1

   # Elais
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),CIRCLE('ICRS',242.5,54.0,2.5))=1

   # Aegis
   # SELECT * FROM public.igsl_source  WHERE CONTAINS(POINT('ICRS',public.igsl_source.ra,public.igsl_source.dec),CIRCLE('ICRS',214.50,52.50,2.0))=1


   $STILTS  tcatn nin=12 in1=$DATADIR/IGSL/Fall1.fits in2=$DATADIR/IGSL/Fall2a.fits in3=$DATADIR/IGSL/Fall2b.fits in4=$DATADIR/IGSL/Fall3.fits \
      in5=$DATADIR/IGSL/Spring1.fits in6=$DATADIR/IGSL/Spring2.fits in7=$DATADIR/IGSL/Spring3.fits in8=$DATADIR/IGSL/Spring4.fits in9=$DATADIR/IGSL/Spring5.fits \
      in10=$DATADIR/IGSL/Northern.fits in11=$DATADIR/IGSL/Elais.fits  in12=$DATADIR/IGSL/Aegis.fits  \
      out=$DATADIR/IGSL/IGSL_STARS_HSC_footprint.fits \
      ocmd='select "(source_classification != 2 && classification == false) || (source_classification == 2 && classification == true)"'

   return
}


function getPickles
{
   # http://vizier.u-strasbg.fr/viz-bin/VizieR VI135/15 (all columns)

   $STILTS tcatn nin=7 icmd6='delcols _r' icmd7='delcols _r' in1=$DATADIR/Pickles2010/Fall_1.fits in2=$DATADIR/Pickles2010/Fall_2.fits in3=$DATADIR/Pickles2010/Fall_3.fits in4=$DATADIR/Pickles2010/Spring.fits in5=$DATADIR/Pickles2010/Northern.fits in6=$DATADIR/Pickles2010/Elais.fits in7=$DATADIR/Pickles2010/Aegis.fits \
      out=$DATADIR/Pickles2010/Pickles2010_HSC_footprint.fits ocmd='addcol source_id $0; addcol G_Gaia "gfmag-0.0940-0.5310*(gfmag-ifmag)-0.0974*pow(gfmag-ifmag,2.0) + 0.0052*pow(gfmag-ifmag, 3.0)"'

   $STILTS tpipe  in=$DATADIR/Pickles2010/pipeTest.fits \
      out=$DATADIR/Pickles2010/Pickles2010_pipeTest.fits \
      cmd='addcol source_id $0; addcol G_Gaia "gfmag-0.0940-0.5310*(gfmag-ifmag)-0.0974*pow(gfmag-ifmag,2.0) + 0.0052*pow(gfmag-ifmag, 3.0)"'

   return
}


function getHSC
{
   # python scripts/hscSspQuery.py config/tract.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/tract_9376.fits
   # python scripts/hscSspQuery.py config/tract_randoms.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/tract_randoms_9376.fits

   # python scripts/hscSspQuery.py config/wide_randoms.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/wide_randoms.fits
   # python scripts/hscSspQuery.py config/deep_randoms.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/deep_randoms.fits
   # python scripts/hscSspQuery.py config/udeep_randoms.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/udeep_randoms.fits

   # python scripts/hscSspQuery.py config/wide_saturated.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/wide_bright.fits
   # python scripts/hscSspQuery.py config/deep_saturated.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/deep_bright.fits
   # python scripts/hscSspQuery.py config/udeep_saturated.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/udeep_bright.fits

   # $STILTS tcatn nin=3  in1=$DATADIR/HSC-SSP/wide_randoms.fits  in2=$DATADIR/HSC-SSP/deep_randoms.fits  in3=$DATADIR/HSC-SSP/udeep_randoms.fits out=$DATADIR/HSC-SSP/all_randoms.fits
   # $STILTS tcatn nin=3  in1=$DATADIR/HSC-SSP/wide_bright.fits  in2=$DATADIR/HSC-SSP/deep_bright.fits  in3=$DATADIR/HSC-SSP/udeep_bright.fits out=$DATADIR/HSC-SSP/all_bright.fits

   python scripts/hscSspQuery.py config/wide_i22_5.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/wide_i22_5.fits
   # python scripts/hscSspQuery.py config/wide_i22_5_randoms.sql -r dr1 -u "couponj" -D -M  --format fits > $DATADIR/HSC-SSP/wide_i22_5_randoms.fits

   return
}

function getCFHTLenS
{

   $STILTS tcatn nin=4 \
      in1=$HOME/data/CFHTLenS/randoms/W1_ran_MASK_all.ascii ifmt1=ascii icmd1='colmeta -name ra col1; colmeta -name dec col2; colmeta -name MASK col3; select "(MASK == 0 || MASK == 1)"' \
      in2=$HOME/data/CFHTLenS/randoms/W2_ran_MASK_all.ascii ifmt2=ascii icmd2='colmeta -name ra col1; colmeta -name dec col2; colmeta -name MASK col3; select "(MASK == 0 || MASK == 1)"' \
      in3=$HOME/data/CFHTLenS/randoms/W3_ran_MASK_all.ascii ifmt3=ascii icmd3='colmeta -name ra col1; colmeta -name dec col2; colmeta -name MASK col3; select "(MASK == 0 || MASK == 1)"' \
      in4=$HOME/data/CFHTLenS/randoms/W4_ran_MASK_all.ascii ifmt4=ascii icmd4='colmeta -name ra col1; colmeta -name dec col2; colmeta -name MASK col3; select "(MASK == 0 || MASK == 1)"' \
      out=$DATADIR/CFHTLenS/CFHTLenS_i22_5_randoms.fits

   $STILTS tcatn nin=4 \
      in1=$HOME/data/CFHTLenS/release/W1_release.fits icmd1='colmeta -name ra ALPHA_J2000; colmeta -name dec DELTA_J2000; addcol MAGERR_yd toDouble(MAGERR_y); delcols MAGERR_y; addcol MAG_iy "MAG_i > 0.0 ? MAG_i : MAG_y"; select "10.0 < MAG_iy && MAG_iy < 22.5 && (MASK == 0 || MASK == 1)"' \
      in2=$HOME/data/CFHTLenS/release/W2_release.fits icmd2='colmeta -name ra ALPHA_J2000; colmeta -name dec DELTA_J2000; addcol MAGERR_yd toDouble(MAGERR_y); delcols MAGERR_y; addcol MAG_iy "MAG_i > 0.0 ? MAG_i : MAG_y"; select "10.0 < MAG_iy && MAG_iy < 22.5 && (MASK == 0 || MASK == 1)"' \
      in3=$HOME/data/CFHTLenS/release/W3_release.fits icmd3='colmeta -name ra ALPHA_J2000; colmeta -name dec DELTA_J2000; addcol MAGERR_yd toDouble(MAGERR_y); delcols MAGERR_y; addcol MAG_iy "MAG_i > 0.0 ? MAG_i : MAG_y"; select "10.0 < MAG_iy && MAG_iy < 22.5 && (MASK == 0 || MASK == 1)"' \
      in4=$HOME/data/CFHTLenS/release/W4_release.fits icmd4='colmeta -name ra ALPHA_J2000; colmeta -name dec DELTA_J2000; addcol MAGERR_yd toDouble(MAGERR_y); delcols MAGERR_y; addcol MAG_iy "MAG_i > 0.0 ? MAG_i : MAG_y"; select "10.0 < MAG_iy && MAG_iy < 22.5 && (MASK == 0 || MASK == 1)"' \
      out=$DATADIR/CFHTLenS/CFHTLenS_i22_5.fits

   return
}



function getPanSTARRS
{

   # queries launched from http://mastweb.stsci.edu/ps1casjobs/SubmitJob.aspx

   # decimal degrees
   # fall1 dbo.fGetObjFromRectEq(329.00,-02.00,360.0,08.00)
   # fall2 dbo.fGetObjFromRectEq(0.0,-02.00,41.0,08.00)
   # fall3 dbo.fGetObjFromRectEq(26.50,-08.00,41.00,-02.00)
   # spring1 dbo.fGetObjFromRectEq(126.5,-03.00,180.0,06.00)
   # spring2 dbo.fGetObjFromRectEq(180,-03.00,226.0,06.00)
   # Northern dbo.fGetObjFromRectEq(199.0,41.5,251.00,45.00)
   # Elais dbo.fGetNearbyObjEq(242.5,55.0,150.0)
   # Aegis dbo.fGetNearbyObjEq(214.5,52.5,120.0)

   # use query below with these changes:
   # INTO mydb.field_name
   # AND o.gMeanPSFMag < 20.0 OR o.rMeanPSFMag  < 20.0 OR o.iMeanPSFMag < 20.0

   # $STILTS tcatn nin=8 in1=$DATADIR/SDSS/fall1.fits in2=$DATADIR/SDSS/fall2.fits in3=$DATADIR/SDSS/fall3.fits in4=$DATADIR/SDSS/spring1.fits in5=$DATADIR/SDSS/spring2.fits in6=$DATADIR/SDSS/northern.fits in7=$DATADIR/SDSS/elais.fits in8=$DATADIR/SDSS/aegis.fits out=$DATADIR/SDSS/SDSS_HSC_footprint.fits

   # return

   # see https://confluence.stsci.edu/display/PANSTARRS/PS1+Sample+queries

   QUERY="
   SELECT  o.objID,
      ot.raStack, ot.decStack, ot.raMean, ot.decMean,
      ot.ng,  o.gMeanPSFMag,o.gMeanPSFMagErr,o.gMeanKronMag,o.gMeanKronMagErr,
      ot.nr,  o.rMeanPSFMag,o.rMeanPSFMagErr,o.rMeanKronMag,o.rMeanKronMagErr,
      ot.ni,  o.iMeanPSFMag,o.iMeanPSFMagErr,o.iMeanKronMag,o.iMeanKronMagErr,
      ot.nz,  o.zMeanPSFMag,o.zMeanPSFMagErr,o.zMeanKronMag,o.zMeanKronMagErr,
      ot.ny,  o.yMeanPSFMag,o.yMeanPSFMagErr,o.yMeanKronMag,o.yMeanKronMagErr,
      o.gQfPerfect,o.rQfPerfect,o.iQfPerfect,o.zQfPerfect,o.yQfPerfect,
      ot.qualityFlag,ot.objInfoFlag,
      soa.gpsfLikelihood, soa.gKronRad, soa.rKronRad, soa.iKronRad,
      sov.ginfoFlag,sov.rinfoFlag,sov.iinfoFlag,sov.zinfoFlag,sov.yinfoFlag

      INTO mydb.tract_9376
      FROM MeanObject AS o
      JOIN dbo.fGetObjFromRectEq(221.48148, -1.4876, 222.96296, -5.3e-7) r ON r.objid = o.objID
      JOIN ObjectThin AS ot ON ot.objID = o.objID
      LEFT JOIN StackObjectAttributes AS soa ON soa.objID = o.objID
      LEFT JOIN StackObjectView AS sov ON sov.objID = o.objID

      WHERE ot.ni >= 3
      AND ot.ng >= 3
      AND ot.nr >= 3
   "

   cd config/PanSTARRS

	$CASJOB run -t "PanSTARRS_dr1" "$QUERY"
   $CASJOB extract -b mydb.tract_9376 -F -type fits # -d $DATADIR/PanSTARRS/ (NOT WORKING download directly from http://mastweb.stsci.edu/ps1casjobs/Output.aspx)

   cd -

   return
}


function getSDSS
{

   # queries launched from https://skyserver.sdss.org/CasJobs/SubmitJob.aspx
   # downloaded from http://www.voservices.net/skyquery/Apps/MyDb/Tables.aspx?dataset=MYDB

   # decimal degrees
   # fall1 dbo.fGetObjFromRectEq(329.00,-02.00,360.0,08.00)
   # fall2 dbo.fGetObjFromRectEq(0.0,-02.00,41.0,08.00)
   # fall3 dbo.fGetObjFromRectEq(26.50,-08.00,41.00,-02.00)
   # Spring1 dbo.fGetObjFromRectEq(126.5,-03.00,180.0,06.00)
   # Spring2 dbo.fGetObjFromRectEq(180,-03.00,226.0,06.00)
   # Northern dbo.fGetObjFromRectEq(199.0,41.5,251.00,45.00)
   # Elais dbo.fGetNearbyObjEq(242.5,55.0,150.0)
   # Aegis dbo.fGetNearbyObjEq(214.5,52.5,120.0)
   # pipeTest dbo.fGetObjFromRectEq(317.0,-3.0, 324.0, 4.0)

   # INTO mydb.field_name
   # where obj.psfMag_g < 20.0 OR obj.psfMag_r < 20.0 OR obj.psfMag_i < 20.
   # $STILTS tcatn nin=8 in1=$DATADIR/SDSS/fall1.fits in2=$DATADIR/SDSS/fall2.fits in3=$DATADIR/SDSS/fall3.fits in4=$DATADIR/SDSS/spring1.fits in5=$DATADIR/SDSS/spring2.fits in6=$DATADIR/SDSS/northern.fits in7=$DATADIR/SDSS/elais.fits in8=$DATADIR/SDSS/aegis.fits out=$DATADIR/SDSS/SDSS_HSC_footprint.fits

   # cp $DATADIR/SDSS/pipeTest.fits $DATADIR/SDSS/SDSS_pipeTest.fits

   # tract 9376
   if false; then

      QUERY="
      SELECT obj.objID, obj.ra,obj.dec,
        obj.type, obj.clean, obj.probPSF,
        obj.psfMag_u, obj.psfMag_g, obj.psfMag_r, obj.psfMag_i, obj.psfMag_z,
        obj.psfMagErr_u, obj.psfMagErr_g, obj.psfMagErr_r, obj.psfMagErr_i, obj.psfMagErr_z,
        obj.cModelMag_u, obj.cModelMag_g, obj.cModelMag_r, obj.cModelMag_i, obj.cModelMag_z,
        obj.cModelMagErr_u, obj.cModelMagErr_g, obj.cModelMagErr_r, obj.cModelMagErr_i, obj.cModelMagErr_z,
        obj.extinction_u, obj.extinction_g, obj.extinction_r, obj.extinction_i, obj.extinction_z,
        (obj.flags_u & dbo.fPhotoFlags('SATUR_CENTER')) as saturCenter_u,
        (obj.flags_g & dbo.fPhotoFlags('SATUR_CENTER')) as saturCenter_g,
        (obj.flags_r & dbo.fPhotoFlags('SATUR_CENTER')) as saturCenter_r,
        (obj.flags_i & dbo.fPhotoFlags('SATUR_CENTER')) as saturCenter_i,
        (obj.flags_z & dbo.fPhotoFlags('SATUR_CENTER')) as saturCenter_z,

        (obj.flags_u & dbo.fPhotoFlags('INTERP_CENTER')) as interpCenter_u,
        (obj.flags_g & dbo.fPhotoFlags('INTERP_CENTER')) as interpCenter_g,
        (obj.flags_r & dbo.fPhotoFlags('INTERP_CENTER')) as interpCenter_r,
        (obj.flags_i & dbo.fPhotoFlags('INTERP_CENTER')) as interpCenter_i,
        (obj.flags_z & dbo.fPhotoFlags('INTERP_CENTER')) as interpCenter_z

        FROM PhotoPrimary as obj
        JOIN dbo.fGetObjFromRectEq(221.48148, -1.4876, 222.96296, -5.3e-7) r ON r.objid = obj.objID
        INTO mydb.tract_9376
      "

     cd config/SDSS

     $CASJOB run -t "dr13/1" "$QUERY"
     $CASJOB extract -b mydb.tract_9376 -F -type fits

     cd -

  fi

  return

}



# ----------------------------------------------------------------------------------- #
# utils
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------------- #

main $@
