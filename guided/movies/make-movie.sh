#!/bin/bash

mogrify -format jpg *.png
ffmpeg -r 5 -i HealAstr%05d.jpg movie.mov
rm *.jpg
