#! /bin/bash
set -o nounset

SOURCE=text.tex
#EXT=pdf

main() {



  DIR=$( echo $HOME/projects/HSC/masks/brightObjectMasks/plots $HOME/data/HSC/brightObjectMasks/crossCorr/mag_{03.50,04.50,05.50,06.50,07.50,09.50,11.50,13.50,15.50,17.50}_images )
  FIGURES=$( grep includegraphics $SOURCE | awk -F"{{" '{print substr($NF,1,length($NF)-1)}')
  for f in $FIGURES; do
    for d in $DIR; do

#      if [ -e $d/${f}.$EXT ]; then
#        cp $d/${f}.$EXT figures/
#      fi

      FILE=$d/${f/\}}
      if [  -e $FILE ]; then
        echo "Importing $FILE..."
        cp $FILE figures/
      fi

    done
  done


}

# ----------------------------------------------------------------------------------- #
# functions
# ----------------------------------------------------------------------------------- #

function sayHello
{

  echo "Hello"

}


# ----------------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------------- #

main
