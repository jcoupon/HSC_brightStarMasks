#! /bin/bash
set -o nounset

# ----------------------------------------------------------------------------------- #
# people contacted for checks:
#     - Masao Hayashi
#     - Alexie Leauthaud
#     - Surhud More
#     - Masayuki Akiyama
#     - Hisakazu Uchiyama
#     - Elinor Medezinski
#     - Song Huang
#     - Masamune Oguri
#     - Ryoma Murata
#     - Yutaka Komiyama
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# options
# ----------------------------------------------------------------------------------- #

source scripts/config.sh

main() {

   # ----------------------------------------------------------------------------------- #
   # get data
   # ----------------------------------------------------------------------------------- #

   # scripts/getdata.sh

   # ----------------------------------------------------------------------------------- #
   # build star catalogue and perform first tests
   # ----------------------------------------------------------------------------------- #

   # assembleCats
   # doPlots
   # checkPreviousCat

   # ----------------------------------------------------------------------------------- #
   # build mask
   # ----------------------------------------------------------------------------------- #

   MAG=(   03.50  04.50  05.50 06.50 07.50 09.50 11.50 13.50 15.50 17.50 )
   DMAG=(  1.0    1.0    1.0   0.5   0.5   0.1   0.04  0.01  0.004 0.002 )
   LIMIT=( 1200.0 1200.0 800.0 400.0 400.0 250.0 150.0 100.0 60.0  40.0 )

   # makeCzakonMask
   # getInfo
   # crossMatch
   # makeMask
   # getImages
   # makeRegionFiles
   # doTPCFTests
   # plotCAMIRA
   # plotPSFwGhost
   # plotExample
   # makeRegionFilesTest

   check_tracts

   # for version name see http://www.astro.wisc.edu/~dolan/constellations/extra/brightest.html
   # makeRelease Arcturus

   return
}

# ----------------------------------------------------------------------------------- #
# functions
# ----------------------------------------------------------------------------------- #


function check_tracts
{

  scripts/buildMasks.py get_tract_coord -skyMap info/skyMapHSC-SSP.pickle \
    -i info/tracts_s18a.lis.txt -o info/tracts_s18a.lis.coord.fits

  return

  scripts/buildMasks.py get_tract_coord -skyMap info/skyMapHSC-SSP.pickle \
    -i info/tracts_arcturus_20180226.lis.txt \
    -o info/tracts_arcturus_20180226.lis.coord.fits


  return
}



function assembleCats
# Assemble Gaia and Tycho-2. When the source is duplicated,
# keep Gaia position and magnitude.
# Also build a matched Gaia+Tycho-2 for tests
{

   # for pipeline tests
   if false; then


      # assemble Gaia and Tycho-2
      $STILTS tmatch2 join=1and2 in1=$DATADIR/Pickles2010/Pickles2010_pipeTest.fits in2=$DATADIR/Gaia/Gaia/Gaia_pipeTest.fits out=$DATADIR/Pickles2010/Pickles2010_pipeTest_Gaia.fits \
         matcher=sky params=1 values1='RAJ2000 DEJ2000' values2='ra dec' suffix1="" suffix2="_Gaia"

      $STILTS tmatch2 join=1not2 in1=$DATADIR/Pickles2010/Pickles2010_pipeTest.fits in2=$DATADIR/Gaia/Gaia_pipeTest.fits out=$DATADIR/Pickles2010/Pickles2010_pipeTest_no_Gaia_match.fits \
         matcher=sky params=1 values1='RAJ2000 DEJ2000' values2='ra dec'

      $STILTS tcatn nin=2 in1=$DATADIR/Gaia/Gaia_pipeTest.fits in2=$DATADIR/Pickles2010/Pickles2010_pipeTest_no_Gaia_match.fits \
         icmd1='addcol origin \"Gaia\";select "duplicated_source != true"; colmeta -name G_Gaia phot_g_mean_mag; keepcols "source_id ra dec G_Gaia origin"' icmd2='addcol origin \"Tycho2\"; colmeta -name ra RAJ2000; colmeta -name dec DEJ2000; keepcols "source_id ra dec G_Gaia origin"' \
         out=$DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest.fits

      # HSC wide footprint
      $STILTS tmatchn matcher=skyerr params=1.0 join1=always nin=2 \
         in1=$DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest.fits values1='ra dec 1.0' suffix1="" \
         in2=$DATADIR/SDSS/SDSS_pipeTest.fits  values2='ra dec 1.0' suffix2="_SDSS" icmd2='select "clean == 1"; addcol G_Gaia "psfMag_g-0.0940-0.5310*(psfMag_g-psfMag_i)-0.0974*pow(psfMag_g-psfMag_i,2.0) + 0.0052*pow(psfMag_g-psfMag_i, 3.0)"; addcol extended_g "0.1 < psfMag_g-cModelMag_g && psfMag_g-cModelMag_g < 20.0 ? 1:0"; addcol extended_r "0.1 < psfMag_r-cModelMag_r && psfMag_r-cModelMag_r < 20.0 ? 1:0"; addcol extended_i "0.1 < psfMag_i-cModelMag_i && psfMag_i-cModelMag_i < 20.0 ? 1:0"; addcol extended "(saturCenter_g == 0 && saturCenter_r == 0 && saturCenter_i == 0) && (extended_g == 1 || extended_r == 1 || extended_i == 1) ? 1:0"; keepcols "ra dec G_Gaia extended"' \
         out=$DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest_SDSS.fits ocmd='delcols "ra_SDSS dec_SDSS"'

      # this is the one catalogue
      $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest_SDSS.fits out=$DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest_SDSS_stars_pure.fits cmd='select "G_Gaia < 18.0 && !(G_Gaia > 14 && extended == 1)"; addcol r_Czakon_deg  "(200.0*pow(10.0, 0.25*(7.0 - G_Gaia)) + 12.0*pow(10.0, 0.05*(16.0 - G_Gaia))) / 3600.0"; delcols extended'
      exit
   fi


   if false; then
      # assemble Gaia and Tycho-2
      $STILTS tmatch2 join=1and2 in1=$DATADIR/Pickles2010/Pickles2010_HSC_footprint.fits in2=$DATADIR/Gaia/Gaia_HSC_footprint.fits out=$DATADIR/Pickles2010/Pickles2010_HSC_footprint_Gaia.fits \
         matcher=sky params=1 values1='RAJ2000 DEJ2000' values2='ra dec' suffix1="" suffix2="_Gaia"

      $STILTS tmatch2 join=1not2 in1=$DATADIR/Pickles2010/Pickles2010_HSC_footprint.fits in2=$DATADIR/Gaia/Gaia_HSC_footprint.fits out=$DATADIR/Pickles2010/Pickles2010_HSC_footprint_no_Gaia_match.fits \
         matcher=sky params=1 values1='RAJ2000 DEJ2000' values2='ra dec'

      $STILTS tcatn nin=2 in1=$DATADIR/Gaia/Gaia_HSC_footprint.fits in2=$DATADIR/Pickles2010/Pickles2010_HSC_footprint_no_Gaia_match.fits \
         icmd1='addcol origin \"Gaia\";select "duplicated_source != true"; colmeta -name G_Gaia phot_g_mean_mag; keepcols "source_id ra dec G_Gaia origin"' icmd2='addcol origin \"Tycho2\"; colmeta -name ra RAJ2000; colmeta -name dec DEJ2000; keepcols "source_id ra dec G_Gaia origin"' \
         out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint.fits

      # HSC wide footprint
      $STILTS tmatchn matcher=skyerr params=1.0 join1=always nin=4 \
         in1=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint.fits values1='ra dec 1.0' suffix1="" \
         in2=$DATADIR/SDSS/SDSS_HSC_footprint.fits  values2='ra dec 1.0' suffix2="_SDSS" icmd2='select "clean == 1"; addcol G_Gaia "psfMag_g-0.0940-0.5310*(psfMag_g-psfMag_i)-0.0974*pow(psfMag_g-psfMag_i,2.0) + 0.0052*pow(psfMag_g-psfMag_i, 3.0)"; addcol extended_g "0.1 < psfMag_g-cModelMag_g && psfMag_g-cModelMag_g < 20.0 ? 1:0"; addcol extended_r "0.1 < psfMag_r-cModelMag_r && psfMag_r-cModelMag_r < 20.0 ? 1:0"; addcol extended_i "0.1 < psfMag_i-cModelMag_i && psfMag_i-cModelMag_i < 20.0 ? 1:0"; addcol extended "(saturCenter_g == 0 && saturCenter_r == 0 && saturCenter_i == 0) && (extended_g == 1 || extended_r == 1 || extended_i == 1) ? 1:0"; keepcols "ra dec G_Gaia extended"' \
         in3=$DATADIR/HSC-SSP/all_bright.fits values3='ra dec 1.0' suffix3="_HSC"  icmd3='select "parent_id == 0"; addcol G_Gaia "gmag_psf-0.0940-0.5310*(gmag_psf-imag_psf)-0.0974*pow(gmag_psf-imag_psf,2.0) + 0.0052*pow(gmag_psf-imag_psf, 3.0)"; addcol extended "!iflags_pixel_saturated_center && 0.01 < imag_psf - icmodel_mag && imag_psf - icmodel_mag  < 20.0 ? 1:0 "; keepcols "ra dec G_Gaia extended"' \
         in4=$DATADIR/HSC-SSP/all_randoms.fits values4='ra dec 120.0' suffix4="_ran"  icmd4='addcol ipsf_FWHM "ipsf_pix*0.17*2.0*sqrt(2.0)"; keepcols "ra dec ipsf_FWHM"' \
         out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_S16A_SDSS.fits ocmd='delcols "ra_HSC dec_HSC ra_SDSS dec_SDSS ra_ran dec_ran "'

   fi

   # this is the one catalogue
   $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_S16A_SDSS.fits out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits cmd='select "G_Gaia < 18.0 && !(G_Gaia > 14 && extended_SDSS == 1)"; addcol r_Czakon_deg  "(200.0*pow(10.0, 0.25*(7.0 - G_Gaia)) + 12.0*pow(10.0, 0.05*(16.0 - G_Gaia))) / 3600.0"; delcols extended_SDSS'

   # 1.47082902202 % extended between 14 and 18

   # objects extended in HSC
   $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure_HSC_extended.csv ofmt=csv cmd='select "extended_HSC == 1"'

   # tract 9376 - pure
   $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure_tract_9376.fits cmd='select "221.48148 < ra && ra < 222.96296 && -1.4876 < dec && dec < -5.3e-7"'

   # tract 9376 - all
   $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint.fits out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_tract_9376.fits cmd='select "221.48148 < ra && ra < 222.96296 && -1.4876 < dec && dec < -5.3e-7"'

   $STILTS tmatchn matcher=sky params=1.0 join1=always nin=4 \
      in1=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_tract_9376.fits values1='ra dec' suffix1="" \
      in2=$DATADIR/HSC-SSP/tract_9376.fits values2='ra dec' suffix2="_HSC"  icmd2='select "parent_id == 0"; addcol G_Gaia "gmag_psf-0.0940-0.5310*(gmag_psf-imag_psf)-0.0974*pow(gmag_psf-imag_psf,2.0) + 0.0052*pow(gmag_psf-imag_psf, 3.0)"; addcol extended "!iflags_pixel_saturated_center && 0.01 < imag_psf - icmodel_mag && imag_psf - icmodel_mag  < 20.0 ? 1:0 "; keepcols "ra dec G_Gaia extended"' \
      in3=$DATADIR/SDSS/tract_9376.fits  values3='ra dec' suffix3="_SDSS" icmd3='select "clean == 1"; addcol G_Gaia "psfMag_g-0.0940-0.5310*(psfMag_g-psfMag_i)-0.0974*pow(psfMag_g-psfMag_i,2.0) + 0.0052*pow(psfMag_g-psfMag_i, 3.0)"; addcol extended_g "0.1 < psfMag_g-cModelMag_g && psfMag_g-cModelMag_g < 20.0 ? 1:0"; addcol extended_r "0.1 < psfMag_r-cModelMag_r && psfMag_r-cModelMag_r < 20.0 ? 1:0"; addcol extended_i "0.1 < psfMag_i-cModelMag_i && psfMag_i-cModelMag_i < 20.0 ? 1:0"; addcol extended "(saturCenter_g == 0 && saturCenter_r == 0 && saturCenter_i == 0) && (extended_g == 1 || extended_r == 1 || extended_i == 1) ? 1:0"; keepcols "ra dec G_Gaia extended"' \
      in4=$DATADIR/PanSTARRS/PanSTARRS_HSC_tract_9376.fits values4='raMean decMean' suffix4="_PS1" icmd4='addcol extended "0.0 < iMeanPSFMag-iMeanKronMag && iMeanPSFMag-iMeanKronMag < 20.0 ? 1:0";  addcol G_Gaia "gMeanPSFMag-0.0940-0.5310*(gMeanPSFMag-iMeanPSFMag)-0.0974*pow(gMeanPSFMag-iMeanPSFMag,2.0) + 0.0052*pow(gMeanPSFMag-iMeanPSFMag, 3.0)"; keepcols "raMean decMean G_Gaia extended"' \
      out=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_tract_9376_HSC-SSP_SDSS_PS1.fits ocmd='delcols "ra_HSC dec_HSC ra_SDSS dec_SDSS raMean decMean "'

   return
}

function doPlots
{


   # $VENICE -m $DATADIR/reg/HSC-SSP_brightStarMask.reg -r -coord spher -f inside -cd -npart 80000 -o "!$DATADIR/reg/inside.fits"

   # scripts/plotResults.py skyMap -i $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits -o plots/starDensity.pdf -l density
   # scripts/plotResults.py skyMap -i $DATADIR/reg/inside.fits -o plots/maskedFraction.pdf -l masked


   if false; then
      scripts/plotResults.py plotSeeingDist -i $DATADIR/HSC-SSP/wide_randoms.fits -o plots/seeingDist_wide.pdf -t "Wide"
      scripts/plotResults.py plotSeeingDist -i $DATADIR/HSC-SSP/deep_randoms.fits -o plots/seeingDist_deep.pdf -t "Deep" --nolegend
      scripts/plotResults.py plotSeeingDist -i $DATADIR/HSC-SSP/udeep_randoms.fits -o plots/seeingDist_udeep.pdf -t "Ultra Deep" --nolegend

      scripts/plotResults.py plotSaturVsSeeing -i $DATADIR/HSC-SSP/wide_bright.fits  -o plots/saturVsSeeing_wide.pdf -s "(iflags_pixel_saturated_center == True) & (parent_id == 0)" -t "Wide" --nolegend
      scripts/plotResults.py plotSaturVsSeeing -i $DATADIR/HSC-SSP/deep_bright.fits  -o plots/saturVsSeeing_deep.pdf -s "(iflags_pixel_saturated_center == True) & (parent_id == 0)" -t "Deep" --nolegend
      scripts/plotResults.py plotSaturVsSeeing -i $DATADIR/HSC-SSP/udeep_bright.fits  -o plots/saturVsSeeing_udeep.pdf -s "(iflags_pixel_saturated_center == True) & (parent_id == 0)" -t "Ultra Deep" --nolegend
   fi


   if false; then
      stiff "$DATADIR/crossCorr/mag_03.50_images/HSC-R/13-cutout-HSC-R-9805-s16a_wide.fits[1]" "$DATADIR/crossCorr/mag_03.50_images/HSC-I/13-cutout-HSC-I-9805-s16a_wide.fits[1]" "$DATADIR/crossCorr/mag_03.50_images/HSC-Z/13-cutout-HSC-Z-9805-s16a_wide.fits[1]" \
         -NEGATIVE N -BINNING 10,10  -OUTFILE_NAME plots/saturated_star.tif # -COMPRESSION_TYPE JPEG

      stiff "$DATADIR/crossCorr/mag_03.50_images/HSC-R/13-cutout-HSC-R-9805-s16a_wide.fits[1]" "$DATADIR/crossCorr/mag_03.50_images/HSC-I/13-cutout-HSC-I-9805-s16a_wide.fits[1]" "$DATADIR/crossCorr/mag_03.50_images/HSC-Z/13-cutout-HSC-Z-9805-s16a_wide.fits[1]" \
         -NEGATIVE Y -BINNING 10,10  -OUTFILE_NAME plots/saturated_star_neg.tif # -COMPRESSION_TYPE JPEG
   fi

   # scripts/plotResults.py plotMagDist -i $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint.fits,$DATADIR/Pickles2010/Pickles2010_HSC_footprint.fits -o plots/magDist_GaiaTycho-2.pdf -k "G_Gaia","G_Gaia" -s "origin == 'Gaia'","" -l "Gaia","Tycho-2" -t "HSC footprint"
   # scripts/plotResults.py plotMagDiff -i $DATADIR/Pickles2010/Pickles2010_HSC_footprint_Gaia.fits -o plots/magDiff_Tycho-2_emulated.pdf  -t "Tycho-2 emulated versus real Gaia magnitude"

   if true; then
      scripts/plotResults.py plotGalFrac -i $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_S16A_SDSS.fits,$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_S16A_SDSS.fits -s "G_Gaia_HSC > 0.0","G_Gaia_SDSS > 0.0" -o plots/galfrac.pdf
      # scripts/plotResults.py plotMagDist -i $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_SDSS_PS1.fits,$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_SDSS_PS1.fits,$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_HSC-SSP_SDSS_PS1.fits -o plots/magDist_GalaxyFrac.pdf -k "G_Gaia","G_Gaia","G_Gaia" -s "(extended_HSC==1) & (210.0 < ra) & (ra < 225) & (-2.33 < dec) & (dec < 1.0)","(extended_SDSS==1) & (210.0 < ra) & (ra < 225) & (-2.33 < dec) & (dec < 1.0)","(210.0 < ra) & (ra < 225) & (-2.33 < dec) & (dec < 1.0)" -l "HSC extended","SDSS extended","All"  -t "HSC/SDSS/Gaia common"

   fi

   return
}


function checkPreviousCat
{


   # $STILTS tmatch1 action=keep1 matcher=sky params=0.1 $DATADIR/S16A/all.ascii ifmt=ascii out=$DATADIR/S16A/all_uniq.fits values='ra dec'

   #$STILTS tmatchn nin=3 matcher=sky params=1 join1=always in1=$DATADIR/S16A/all_uniq.fits values1='ra dec' icmd1='select "221.48148<ra&&ra<222.96296 && -1.4876<dec&&dec<-5.3e-7"' \
   #   in2=$DATADIR/SDSS/tract_9376.fits values2='ra dec' suffix2='_SDSS'  icmd2='select "clean == 1 && 221.48148<ra&&ra<222.96296 && -1.4876<dec&&dec<-5.3e-7"; addcol G_Gaia "psfMag_g-0.0940-0.5310*(psfMag_g-psfMag_i)-0.0974*pow(psfMag_g-psfMag_i,2.0) + 0.0052*pow(psfMag_g-psfMag_i, 3.0)"; addcol extended_g "0.1 < psfMag_g-cModelMag_g && psfMag_g-cModelMag_g < 20.0 ? 1:0"; addcol extended_r "0.1 < psfMag_r-cModelMag_r && psfMag_r-cModelMag_r < 20.0 ? 1:0"; addcol extended_i "0.1 < psfMag_i-cModelMag_i && psfMag_i-cModelMag_i < 20.0 ? 1:0"; addcol extended "(saturCenter_g == 0 && saturCenter_r == 0 && saturCenter_i == 0) && (extended_g == 1 || extended_r == 1 || extended_i == 1) ? 1:0"; keepcols "ra dec G_Gaia extended"' \
   #   in3=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits suffix3='' values3='ra dec' icmd3='select "221.48148<ra&&ra<222.96296 && -1.4876<dec&&dec<-5.3e-7"; keepcols "ra dec G_Gaia r_Czakon_arcsec"' \
   #   out=$DATADIR/S16A/all_uniq_tract_9376_SDSSGaiaMatched.fits ocmd='addcol G_Gaia_all "G_Gaia > 0.0 ? G_Gaia : G_Gaia_SDSS" '

   scripts/plotResults.py plotMagDist -i $DATADIR/S16A/all_uniq_tract_9376_SDSSGaiaMatched.fits,$DATADIR/S16A/all_uniq_tract_9376_SDSSGaiaMatched.fits,$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits  -o plots/magDist_S16A.pdf -k "G_Gaia_all","G_Gaia_all","G_Gaia" -s "(G_Gaia>0.0)","(G_Gaia_all>0.0)&(extended==1)","(221.48148<ra) & (ra<222.96296) & (-1.4876<dec) & (dec<-5.3e-7)"  -l "Sirius","Sirius [galaxies]","Arcturus" -t "Tract 9376"

   return
}


function makeCzakonMask
{

   OUTFILE=$DATADIR/CzakonMask/Gaia_Tycho-2_HSC_footprint_stars_pure_Czakon_tract_9376.reg

   echo "# BRIGHT STAR CATALOG: J. Coupon (N. Czakon's radius definition) " > $OUTFILE
   echo "# GENERATED ON: $( date )" >>  $OUTFILE
   echo "# TRACT: None" >> $OUTFILE
   echo "# PATCH: None" >> $OUTFILE
   echo "# FILTER: None" >> $OUTFILE
   echo "wcs; fk5" >> $OUTFILE

   $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits ofmt=ascii cmd='select " 221.48148<ra&&ra<222.96296 && -1.4876<dec&&dec<-5.3e-7"; keepcols "ra dec r_Czakon_deg origin source_id"' \
      | awk '!(/^#/) {printf("circle(%s,%s,%sd) # ID: %s_%s\n", $1, $2, $3, $4, $5)}' >> $OUTFILE



   return
}


function getInfo
{

   ALL=info/crossMatch_all.ascii
   S16A=info/crossMatch_S16A.ascii

   FILEIN=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits

   echo "# mag dmag N MAG_MEAN CZAKON_RADIUS_ARCSEC SEEING_MEAN " > $ALL
   echo "# mag dmag N MAG_MEAN CZAKON_RADIUS_ARCSEC SEEING_MEAN " > $S16A

   for (( m=0; m<${#MAG[@]}; m++ )); do

      # select the star sample
      SELECT='select " '${MAG[m]}'-'${DMAG[m]}'/2.0 < G_Gaia  && G_Gaia < '${MAG[m]}'+'${DMAG[m]}'/2.0 "'

      # record info
      N=$(             $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats | awk -F":" '/Rows/ {print $2}' )
      MAG_MEAN=$(      $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='keepcols "G_Gaia"' | awk -F"|" '/G_Gaia/ {print $3}' )
      R_ARCSEC_MEAN=$( $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='keepcols "r_Czakon_deg"' | awk -F"|" '/r_Czakon_deg/ {print $3*3600.0}' )
      SEEING_MEAN=$(   $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='select "ipsf_FWHM > 0.0"; keepcols "ipsf_FWHM"' | awk -F"|" '/ipsf_FWHM/ {print $3}' )

      echo ${MAG[m]} ${DMAG[m]} $N $MAG_MEAN $R_ARCSEC_MEAN $SEEING_MEAN >> $ALL

      # select the star sample
      SELECT='select " '${MAG[m]}'-'${DMAG[m]}'/2.0 < G_Gaia  && G_Gaia < '${MAG[m]}'+'${DMAG[m]}'/2.0 && !NULL_ipsf_FWHM"'

      # record info
      N=$(             $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats | awk -F":" '/Rows/ {print $2}' )
      MAG_MEAN=$(      $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='keepcols "G_Gaia"' | awk -F"|" '/G_Gaia/ {print $3}' )
      R_ARCSEC_MEAN=$( $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='keepcols "r_Czakon_deg"' | awk -F"|" '/r_Czakon_deg/ {print $3*3600.0}' )
      SEEING_MEAN=$(   $STILTS tpipe $FILEIN cmd="$SELECT" omode=stats cmd='select "ipsf_FWHM > 0.0"; keepcols "ipsf_FWHM"' | awk -F"|" '/ipsf_FWHM/ {print $3}' )

      echo ${MAG[m]} ${DMAG[m]} $N $MAG_MEAN $R_ARCSEC_MEAN $SEEING_MEAN >> $S16A


   done


   return

}

function crossMatch
{

   # naoj password
   USER=couponj
   # read -p "Enter NAOJ password for user $USER: " -s HSC_SSP_CAS_PASSWORD
   export HSC_SSP_CAS_PASSWORD

   for (( m=0; m<${#MAG[@]}; m++ )); do

      # select the star sample
      SELECT='select " '${MAG[m]}'-'${DMAG[m]}'/2.0 < G_Gaia  && G_Gaia < '${MAG[m]}'+'${DMAG[m]}'/2.0 "'
      $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits out=$DATADIR/crossCorr/mag_${MAG[m]}.fits cmd="$SELECT" cmd='select "!NULL_ipsf_FWHM"; keepcols "ra dec G_Gaia r_Czakon_deg ipsf_FWHM"'

      # run query
      # python scripts/hscSspCrossMatch.py $DATADIR/crossCorr/mag_${MAG[m]}.fits --accuracy ${LIMIT[m]} --template config/crossCorr_template.sql > config/crossCorr.sql
      # python scripts/hscSspQuery.py config/crossCorr.sql -r dr1 -u $USER -D -M  --format fits  > $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP.fits

      # run random query -- does not work
      # python scripts/hscSspCrossMatch.py $DATADIR/crossCorr/mag_${MAG[m]}.fits --accuracy $TWO_R_ARCSEC_MEAN --template config/crossCorr_random_template.sql --rerun s16a_wide_random > config/crossCorr_random.sql
      # python scripts/hscSspQuery.py config/crossCorr_random.sql -r dr1 -u $USER -D -M  --format fits  > result_random.fits

      # scripts/stack.py toPixels -i $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP.fits -o $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits

   done

   return
}


function makeMask
{

   OUTFILE=info/maskRadius.ascii
   OUTFILEBESTSEEING=info/maskRadius_bestSeeing.ascii
   OUTFILEWORSTSEEING=info/maskRadius_worstSeeing.ascii

   if false; then
      echo "# mag seeing radiusIsotropic radiusBleedTrail radiusSpike" > $OUTFILE
      echo "# mag seeing radiusIsotropic radiusBleedTrail radiusSpike" > $OUTFILEBESTSEEING
      echo "# mag seeing radiusIsotropic radiusBleedTrail radiusSpike" > $OUTFILEWORSTSEEING
      for (( m=0; m<${#MAG[@]}; m++ )); do

         if true; then

            RESULT=$( scripts/buildMasks.py maskRadius -rmax ${LIMIT[m]}  -mag ${MAG[m]} \
               -i $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits \
               -k "match_distance","match_distance" -s "(parent_id==0)","(deblend_nchild==0)" \
               -o $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_parent.ascii,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_child.ascii )
            echo ${MAG[m]} $RESULT >> $OUTFILE
            echo ${MAG[m]} $RESULT

            RESULT=$( scripts/buildMasks.py maskRadius -rmax ${LIMIT[m]}  -mag ${MAG[m]} -bestSeeing \
               -i $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits \
               -k "match_distance","match_distance" -s "(parent_id==0)","(deblend_nchild==0)" \
               -o $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_parent_bestSeeing.ascii,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_child_bestSeeing.ascii )

            echo ${MAG[m]} $RESULT >> $OUTFILEBESTSEEING
            echo ${MAG[m]} $RESULT "(best seeing)"

         fi


         RESULT=$( scripts/buildMasks.py maskRadius -rmax ${LIMIT[m]}  -mag ${MAG[m]} -worstSeeing \
            -i $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_pix.fits \
            -k "match_distance","match_distance" -s "(parent_id==0)","(deblend_nchild==0)" \
            -o $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_parent_worstSeeing.ascii,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_child_worstSeeing.ascii )

         echo ${MAG[m]} $RESULT >> $OUTFILEWORSTSEEING
         echo ${MAG[m]} $RESULT "(worst seeing)"

      done
   fi

   if false; then
      for (( m=0; m<${#MAG[@]}; m++ )); do
         scripts/plotResults.py plotSourceDensity -o plots/sourceDensityIsotropic_mag_${MAG[m]}.pdf \
            -i $DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_parent.ascii,$DATADIR/crossCorr/mag_${MAG[m]}_HSC-SSP_hist_child.ascii  \
            -l "Parent","Primary" -t "\$G_\mathrm{Gaia}=${MAG[m]}\$" \
            -rmax ${LIMIT[m]} -mag ${MAG[m]} $( [ "$m" != "0" ] && echo "--nolegend" )
      done
   fi

   # scripts/buildMasks.py fitMaskRadius -i $OUTFILE -o ${OUTFILE%.ascii}_bestFit.ascii
   scripts/plotResults.py plotMaskRadius -i $OUTFILE,${OUTFILE%.ascii}_bestFit.ascii,$OUTFILEBESTSEEING,$OUTFILEWORSTSEEING -o plots/maskRadius.pdf

   return
}



function getImages
{

   # naoj password
   USER=couponj
   # read -p "Enter NAOJ password for user $USER: " -s HSC_SSP_CAS_PASSWORD
   export HSC_SSP_CAS_PASSWORD

   RADIUSFILE=info/maskRadius_bestFit.ascii

   for (( m=0; m<${#MAG[@]}; m++ )); do
      for f in HSC-G HSC-R HSC-I HSC-Z HSC-Y; do

   #for m in 6; do
      #  for f in HSC-Y; do

         echo ${MAG[m]} ${f}

         OUTDIR=$DATADIR/crossCorr/mag_${MAG[m]}_images/${f}
         FILELIST=$DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_images_${f}.ascii

         mkdir -p $OUTDIR

         if false; then

            echo "#? filter ra dec sw sh rerun" > $FILELIST
            $STILTS tpipe $DATADIR/crossCorr/mag_${MAG[m]}.fits  ofmt=ascii cmd='keepcols "ra dec"' \
               | awk '!(/^#/) {print "'$f'", $1, $2, "'${LIMIT[m]}'asec", "'${LIMIT[m]}'asec", "s16a_wide " }' >> $FILELIST

            curl https://hscdata.mtk.nao.ac.jp/das_quarry/cgi-bin/quarryImage --form \
               list=@$FILELIST --user $USER:$HSC_SSP_CAS_PASSWORD --insecure  \
               | tar xvf - --strip=1 --directory $OUTDIR

         fi

         R1=$( grep ${MAG[m]} info/maskRadius.ascii | awk '{print $2}' )
         R2=$( grep ${MAG[m]} info/maskRadius.ascii | awk '{print $3}' )
         MAG=$( grep ${MAG[m]} info/maskRadius.ascii | awk '{print $1}' )

         # "\$G_\mathrm{Gaia}=${MAG[m]}\$"

         # scripts/stack.py stackStar -i $OUTDIR -o $DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_${f}.fits
         # scripts/plotResults.py plotStackedStar -i $DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_${f}.fits,$RADIUSFILE  -o $DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_${f}.pdf -t  "\$${MAG[m]}\$" -rmax $R1,$R2,$MAG -l "$f"

         scripts/plotResults.py plotStackedStar -i $DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_${f}.fits,$RADIUSFILE  -o $DATADIR/crossCorr/mag_${MAG[m]}_images/mag_${MAG[m]}_${f}_highRes.pdf -t  "\$${MAG[m]}\$" -rmax $R1,$R2,$MAG -l "$f"


         #return
      done
   done

   return
}



function makeRegionFiles
{

   # tests # ds9 calexp-HSC-Y-9376-5,5.fits -scale mode zscale  -regions load ~/data/HSC/brightStarMasks/tracts/9376/BrightStarMask-9376-5,5-HSC-Y.reg &

   SKYMAP=info/skyMapHSC-SSP.pickle
   RADIUSFILE=info/maskRadius_bestFit.ascii
   STARCAT=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits

   if false; then
      STARCAT=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure_tract_9376.fits
      scripts/buildMasks.py MakeMasterFile -i $STARCAT,$RADIUSFILE -o $DATADIR/reg/HSC-SSP_brightStarMask_tract_9376.reg
      return
   fi

   # scripts/buildMasks.py MakeMasterFile -i $STARCAT,$RADIUSFILE -o $DATADIR/reg/HSC-SSP_brightStarMask.reg
   source $HOME/local/source/hscPipe/env_HSC.sh

   scripts/buildMasks.py makePatchFiles -i $SKYMAP,$STARCAT,$RADIUSFILE -o $DATADIR/brightObjectMasks_reg -basename BrightObjectMask

   cd $DATADIR/brightObjectMasks_reg
   tar czvf tracts.tgz tracts
   tar czvf patches.tgz patches
   cd -


   return

   scripts/buildMasks.py makePatchFiles -i $SKYMAP,$STARCAT,$RADIUSFILE -o $DATADIR/reg

   cd $DATADIR/reg
   tar czvf tracts.tgz tracts
   tar czvf patches.tgz patches
   cd -

   return
}



function makeRegionFilesTest
{
   SKYMAP=info/skyMapTest.pickle
   #SKYMAP=info/skyMapHSC-SSP.pickle
   RADIUSFILE=info/maskRadius_bestFit.ascii
   STARCAT=$DATADIR/starCatalogue/Gaia_Tycho-2_pipeTest_SDSS_stars_pure.fits

   # $STILTS tpipe $DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits out=$STARCAT cmd='replacecol ra "ra-30.0"'

   source $HOME/local/source/hscPipe/env_HSC.sh

   scripts/buildMasks.py makePatchFiles -i $SKYMAP,$STARCAT,$RADIUSFILE -o $DATADIR/brightObjectMasks_pipeTest_reg -basename BrightObjectMask

   return
}



function doTPCFTests
{

   # selection of galaxies and randoms in W4 field
   if false; then

      $VENICE -m $DATADIR/reg/HSC-SSP_brightStarMask.reg -xcol ra -ycol dec -cat $DATADIR/HSC-SSP/wide_i22_5.fits -f all -o "!$DATADIR/HSC-SSP/wide_i22_5_newMask.fits"
      $VENICE -m $DATADIR/reg/HSC-SSP_brightStarMask.reg -xcol ra -ycol dec -cat $DATADIR/HSC-SSP/wide_i22_5_randoms.fits -f all -o "!$DATADIR/HSC-SSP/wide_i22_5_randoms_newMask.fits"

      $STILTS tpipe in=$DATADIR/CFHTLenS/CFHTLenS_i22_5.fits out=$DATADIR/CFHTLenS/CFHTLenS_i22_5_galsW4.fits \
         cmd='select "331.9 < ra && ra < 335.7 && -1.0 < dec && dec < 1.72 && star_flag == 0" '
      $STILTS tpipe in=$DATADIR/CFHTLenS/CFHTLenS_i22_5_randoms.fits out=$DATADIR/CFHTLenS/CFHTLenS_i22_5_randoms_W4.fits \
         cmd='select "331.9 < ra && ra < 335.7 && -1.0 < dec && dec < 1.72"'

      $STILTS tpipe in=$DATADIR/HSC-SSP/wide_i22_5_newMask.fits out=$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits \
         cmd='select "331.9 < ra && ra < 335.7 && -1.0 < dec && dec < 1.72 && iclassification_extendedness != 0"'
      $STILTS tpipe in=$DATADIR/HSC-SSP/wide_i22_5_randoms_newMask.fits out=$DATADIR/HSC-SSP/wide_i22_5_randoms_W4.fits \
         cmd='select "331.9 < ra && ra < 335.7 && -1.0 < dec && dec < 1.72" '

      return
   fi

   # scripts/plotResults.py plotMagDist -i $DATADIR/CFHTLenS/CFHTLenS_i22_5_galsW4.fits,$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits -o plots/magDist_CFHTLenSWide_Arcturus.pdf -k "MAG_i","icmodel_mag-a_i" -s "","flag==1" -l "CFHTLenS","HSC-SSP" -t "Arcturus mask"
   # scripts/plotResults.py plotMagDist -i $DATADIR/CFHTLenS/CFHTLenS_i22_5_galsW4.fits,$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits -o plots/magDist_CFHTLenSWide_Sirius.pdf -k "MAG_i","icmodel_mag-a_i" -s "","(iflags_pixel_bright_object_center==False)" -l "CFHTLenS","HSC-SSP" -t "Sirius mask (old)"

   IMIN=( 17.5 18.5 19.5 20.5 21.5 22.0 )
   IMAX=( 18.5 19.5 20.5 21.5 22.0 22.5 )

   RANGE="-range 0.00005,5.0"

   for (( m=0; m<${#IMIN[@]}; m++ )); do

      echo ${IMIN[m]} ${IMAX[m]}

      scripts/plotResults.py plotWoftheta $( [ "$m" != "0" ] && echo "--nolegend" ) \
         -i $DATADIR/CFHTLenS/CFHTLenS_i${IMIN[m]}_${IMAX[m]}_W4.out,$DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_newMask.out,$DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_oldMask.out,$DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_insideNewMask.out \
         -l "CFHTLenS","Arcturus","Sirius (old)","Inside masks" -o plots/wtheta_i${IMIN[m]}_${IMAX[m]}.pdf -t "\$ ${IMIN[m]} < i < ${IMAX[m]}\$"
      # scripts/plotResults.py plotWoftheta -i $DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_insideNewMask.out,$DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_newMask.out -l "Inside mask","Outside mask" -o plots/wtheta_i${IMIN[m]}_${IMAX[m]}_insideMask.pdf -t "\$ ${IMIN[m]} < i < ${IMAX[m]}\$"
      #exit

      continue

      # CFHTLenS
      $MPIRUN -np $NP $SWOT $RANGE \
         -data1 "$DATADIR/CFHTLenS/CFHTLenS_i22_5_galsW4.fits[${IMIN[m]}<MAG_i&&MAG_i<${IMAX[m]}]" -cols1 ra,dec \
         -ran1 "$DATADIR/CFHTLenS/CFHTLenS_i22_5_randoms_W4.fits" -rancols1 ra,dec -o $DATADIR/CFHTLenS/CFHTLenS_i${IMIN[m]}_${IMAX[m]}_W4.out

      # HSC old mask
      $MPIRUN -np $NP $SWOT $RANGE \
         -data1 "$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits[${IMIN[m]}<icmodel_mag-a_i&&icmodel_mag-a_i<${IMAX[m]}&&.not.iflags_pixel_bright_object_center]" -cols1 ra,dec \
         -ran1 "$DATADIR/HSC-SSP/wide_i22_5_randoms_W4.fits[object_id>-1&&.not.iflags_pixel_bright_object_center]" -rancols1 ra,dec -o $DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_oldMask.out

      # HSC new mask
      $MPIRUN -np $NP $SWOT $RANGE \
         -data1 "$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits[${IMIN[m]}<icmodel_mag-a_i&&icmodel_mag-a_i<${IMAX[m]}&&flag==1]" -cols1 ra,dec \
         -ran1 "$DATADIR/HSC-SSP/wide_i22_5_randoms_W4.fits[flag==1]" -rancols1 ra,dec -o $DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_newMask.out

      # HSC inside new mask
      $MPIRUN -np $NP $SWOT $RANGE \
         -data1 "$DATADIR/HSC-SSP/wide_i22_5_galsW4.fits[${IMIN[m]}<icmodel_mag-a_i&&icmodel_mag-a_i<${IMAX[m]}&&flag==0]" -cols1 ra,dec \
         -ran1 "$DATADIR/HSC-SSP/wide_i22_5_randoms_W4.fits[flag==0]" -rancols1 ra,dec -o $DATADIR/HSC-SSP/wide_i${IMIN[m]}_${IMAX[m]}_W4_insideNewMask.out

   done

   return
}


function plotExample
{

   ds9 -zscale $DATADIR/example/calexp-HSC-G-9799-6,5.fits -region load $DATADIR/example/BrightStarMask-9799-6,5-HSC-R.reg \
      $DATADIR/example/calexp-HSC-R-9799-6,5.fits -region load $DATADIR/example/BrightStarMask-9799-6,5-HSC-R.reg \
      $DATADIR/example/calexp-HSC-I-9799-6,5.fits -region load $DATADIR/example/BrightStarMask-9799-6,5-HSC-R.reg \
      $DATADIR/example/calexp-HSC-Z-9799-6,5.fits -region load $DATADIR/example/BrightStarMask-9799-6,5-HSC-R.reg \
      $DATADIR/example/calexp-HSC-Y-9799-6,5.fits -region load $DATADIR/example/BrightStarMask-9799-6,5-HSC-R.reg &

   return
}

function plotCAMIRA
{
   scripts/plotResults.py CAMIRA \
      -i $DATADIR/CAMIRA/match_sm_arcturus.dat,$DATADIR/CAMIRA/match_sm_sirius.dat -o plots/CAMIRA.pdf \
      -l "Arcturus","Sirius (old)"

   return
}

function plotPSFwGhost
{
   scripts/plotResults.py PSFwGhost -o plots/PSFwGhost.pdf # -t "PSF profile"

   return
}



function makeRelease
{


   RELEASEDIR=$DATADIR/releases/HSC-SSP_brightStarMask_${1}

   rm -rf $RELEASEDIR

   STARCAT=$DATADIR/starCatalogue/Gaia_Tycho-2_HSC_footprint_stars_pure.fits
   REGFILE=$DATADIR/reg/HSC-SSP_brightStarMask.reg

   mkdir -p $RELEASEDIR/stars
   mkdir -p $RELEASEDIR/reg

   cp $REGFILE $RELEASEDIR/reg/masks_all.reg
   cp $DATADIR/reg/tracts.tgz $RELEASEDIR/reg/
   cp $DATADIR/reg/patches.tgz $RELEASEDIR/reg/

   $STILTS tpipe $STARCAT out=$RELEASEDIR/stars/stars_all.fits cmd='delcols "G_Gaia_HSC extended_HSC ipsf_FWHM r_Czakon_deg "'
   cp -r info/venice-?.?.? $RELEASEDIR/

   cp info/README $RELEASEDIR

   cd $DATADIR/releases
   tar czvf HSC-SSP_brightStarMask_${1}.tgz HSC-SSP_brightStarMask_${1}
   cd -

   scp $RELEASEDIR.tgz obs:ftp/brightStarMasks/HSC-SSP

   ssh obs "cd \$HOME/ftp/brightStarMasks/HSC-SSP; ln -sf $( basename $RELEASEDIR.tgz ) HSC-SSP_brightStarMask_latest.tgz  "

   scp info/README obs:ftp/brightStarMasks/HSC-SSP


   return
}

# ----------------------------------------------------------------------------------- #
# utils
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------------- #

main $@
