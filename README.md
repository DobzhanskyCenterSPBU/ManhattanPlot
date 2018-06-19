# ManhattanPlot
Clickable diagram manhatten plot type

The script on Python initating by command line.

All options available on python mp7.py
  -h, --help            show this help message and exit
  --no-log              don't do -log10(p) on the value
  --cols=COLS           zero-based column indexes to get chr, position,
                        p-value respectively e.g. 0,1,2
  --colors=COLORS       cycle through these colors
  --image=IMAGE         save the image to this file. e.g. manhattan7.png
  --title=TITLE         title for the image.
  --ymax=YMAX           max (logged) y-value for plot
  --sep=SEP             data separator, default is [tab]
  --lines               plot the p-values as lines extending from the x-axis
                        rather than points in space. plotting will take longer
                        with this option.
  --dotsactive=DOTSACTIVE
                        Number of clickable top p-value
  --module=MODULE       Module name
  --path=PATH           Path where save the results
  
  Full command like is python mp7.py --cols=0,1,2,3 --colors=rk --dotsactive=19 --path=/home/users/ --module=Modele_NPC myfilecopy.txt

where myfilecopy.txt is a file with records <Chr# Position_insideChr  p-value>.

The file preparing outside the script and connection with the database is progress.

The results of the script are three files:

- graphican png file;

- html file with clicable number of hits (base on -log10(p-value));

- ordered json file with coordinate of hits;

