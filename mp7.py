"""
    %prog [options] files
plot a manhattan plot of the input file(s).
"""

import matplotlib
matplotlib.use('Agg')

import optparse
import sys
from itertools import groupby, cycle
from operator import itemgetter
from matplotlib import pyplot as plt
import numpy as np

def _gen_data(fhs, columns, sep):
    """
    iterate over the files and yield chr, start, pvalue
    """
    for fh in fhs:
        for line in fh:
            if line[0] == "#": continue
            toks = line.split(sep)
            yield toks[columns[0]], int(toks[columns[1]]), float(toks[columns[2]]), str(toks[columns[3]])[0:-1]

def chr_cmp(a, b):
    a = a.lower().replace("_", ""); b = b.lower().replace("_", "")
    achr = a[3:] if a.startswith("chr") else a
    bchr = b[3:] if b.startswith("chr") else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        # X Y
        return cmp(achr, bchr)

def chr_loc_cmp(alocs, blocs):
    return chr_cmp(alocs[0], blocs[0]) or cmp(alocs[1], blocs[1])

def manhattan(fhs, columns, image_path, no_log, colors, sep, title, lines, ymax, dotsactive, module, path):

    xs = []
    ys = []
    cs = []
    colors = cycle(colors)
    xs_by_chr = {}

    last_x = 0
    data = sorted(_gen_data(fhs, columns, sep), cmp=chr_loc_cmp)
#    print "all the data\n", data
#    print "data ends\n"
    
    dataarray = [] 
    dataarray = np.array(data)    
#    print "data as array\n", dataarray
#    print "end of the array\n"

#    print "Again\n"
    for i in range(len(dataarray)):
#	print dataarray[i][0], dataarray[i][1], dataarray[i][2],dataarray[i][3]
	if (float(dataarray[i][2]) < 0):
#		pass
#		print "<0", dataarray[i][2]
		dataarray[i][2] = 0.1
#		print "changed", dataarray[i][2]
#	else:
#		print dataarray[i][0], dataarray[i][1], dataarray[i][2],dataarray[i][3]
 
#    print "end data\n\n\n"
#    print dataarray
   
#    sys.exit(0)
    
#    last_x = 0
#    data = sorted(_gen_data(fhs, columns, sep), cmp=chr_loc_cmp)

    for seqid, rlist in groupby(data, key=itemgetter(0)):
        color = colors.next()
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        xs.extend(region_xs)
        ys.extend([r[2] for r in rlist])

        cs.extend([color] * len(rlist))

        xs_by_chr[seqid] = (region_xs[0] + region_xs[-1]) / 2

        # keep track so that chrs don't overlap.
        last_x = xs[-1]

    xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), cmp=chr_cmp)]

    xs = np.array(xs)

    for i in range(len(ys)):
	if (ys[i] <= 0):
	   ys[i] = 0.1

    ys = np.array(ys) if no_log else -np.log10(ys)


    dataarray = np.array (data)


    for i in range(len(dataarray)):
#    	print dataarray[i][0],dataarray[i][1],dataarray[i][2]
	if (float(dataarray[i][2]) < 0):
		dataarray[i][2] = 0.1

    dataarray = np.array (sorted (dataarray, key=lambda a_entry: float(a_entry[2])  , reverse=False))


#    print "\n\nafter sort\n\n"

#    for i in range(len(dataarray)):
#    	print dataarray[i][0],dataarray[i][1],dataarray[i][2]

    y0p = 0	#46
    y1p = 600	#544
    x1p = 800	#794
    x0p = 0	#90
    xt = 1200000
    yt = 4

#    print (np.amax(xs))
    MAXxs = np.amax(xs)
    MINxs = np.amin(xs)   
#    print ("np.amin(xs)=",MINxs,"np.amax(xs)=",MAXxs)
#    print (np.amax(ys))
    MAXys = np.amax(ys)
    MINys = np.amin(ys)
#    print ("np.amin(ys)=",MINys,"np.amax(ys)=",MAXys)
#    print ("\n")
    
    plt.close()
    pic_width=8
    pic_high=6
    f = plt.figure(facecolor='g',edgecolor='r', figsize=(pic_width,pic_high), dpi=144)
# Define "Canvas" width 8 and high =6 with 144dpi

    ax = f.add_axes ((0.0, 0.0, 1.0, 1.0))  #((0.1, 0.09, 0.88, 0.85))
# Define Area to draw as all the canvas

#    ax.margins(.25)
#    ax.subplots_adjust(0,1,1,0)

#    plt.axes ([1.0, 22.0, 1.0, 10.0])

#    shift = plt.subplots_adjust (left=0.0, right=0.1, top=0.1, bottom=0.0)
#    ylim(0.0,6.0)

    if title is not None:
        plt.title(title)

    ax.set_ylabel('-log10(p-value)')
    ax.set_xlabel('Chr number - in fact, the first field value')

    ax.margins(x=.0, y=0.05)

    if lines:
        ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
    else:
	ax.scatter(xs, ys, s=20, c=cs, alpha=0.8, edgecolors='none')

    N = len(xs)
    htmlfile = path+'mapfile7.html'
#    print htmlfile
    mpfile = open("mapfile7.html","w")

    mpfile.write ("""
<html>
<head>

</head>
<body>

<img src="manhattan7.png" border="0px" width="800px" height="600px" alt="manhattanplot" usemap="#mp" style="margin:0; padding:0;">

""")
    jfile =  open("jasonfile7.txt","w")
    jfile.write("[")

#	for i in range(len(xs)):
    for j in range(int(dotsactive)):
	    i = np.argmax(ys)
	    
#	    plt.plot(xs[i],ys[i],'o',color='yellow')
	    
#	    print "i=",i," ys[i]=",ys[i] ,
#	    mpfile.write("<area shape=\"circle\" coords=\"")

#	    x = (str( -3 + x0p + xs[i] * (x1p - x0p) / MAXxs))
	    x =  (str(6 + x0p + xs[i] * (x1p - x0p) / MAXxs  ))

#	    mpfile.write(x)
#	    mpfile.write(", ")

	    if (ys[i] > 0):
		y =  (str(int( y1p + ys[i] * ((-y1p + y0p) /  (MAXys) ) /1.05 ) + 3 ) ) 
	        x1 = str( x0p + xt * (x1p - x0p) / MAXxs)
	        y1 = str( int (y1p - yt * float(y1p -y0p) / MAXys ))
	    else:
		y = str(0)

#	    yy = str(545 - int(y))
#	    print " y=",y, " yy=",yy
#	    mpfile.write(str(y))
#	    mpfile.write(", ")
#	    mpfile.write("5\"")
#	    mpfile.write(" title=\" ")
#	    mpfile.write(str(i))
#	    mpfile.write(" ")
#	    mpfile.write(data[i][-1])
#	    mpfile.write(" ")
#	    mpfile.write(str(x))
#	    mpfile.write(" ")
#	    mpfile.write(xs[i])
#	    mpfile.write(str(y))
#	    mpfile.write(" ")
#	    mpfile.write(ys[i])
#	    mpfile.write("\" ")
#	    print data[i][-1]
#	    mpfile.write(" href=\" ")
#	    mpfile.write(data[i][-1])
#	    mpfile.write(" ")
#	    mpfile.write(str(x))
#	    mpfile.write(" ")
#	    mpfile.write(str(y))
#	    mpfile.write("\" ")
#	    mpfile.write("  >")

#	    mpfile.write("\n")
	    mpfile.write("<a href=\"\"><img src=\"dot.png\" width=\"5px\" style=\"position:absolute; left:")
	    mpfile.write(str(x))
	    mpfile.write("px; top:")
	    mpfile.write(str(y))
	    mpfile.write("px;\" alt=\"")

	    mpfile.write("module=")
	    mpfile.write(module)
	    mpfile.write(", ")

	    mpfile.write(str(x))
	    mpfile.write(" - ")
	    mpfile.write(str(y))


	    mpfile.write(", ")	    
	    mpfile.write(dataarray[j][0])
	    mpfile.write(", ")
	    mpfile.write(dataarray[j][1])
	    mpfile.write(", ")
	    mpfile.write(dataarray[j][2])

	    mpfile.write("\" ></a>")

#	    mpfile.write("<a href=\"\"><img src=\"dot.png\" width=\"15px\" style=\"position:absolute; left:0; top:0;\"></a>")
#	    mpfile.write("<a href=\"\"><img src=\"dot.png\" width=\"15px\" style=\"position:absolute; left:100px; top:0;\"></a>")
#	    mpfile.write("<a href=\"\"><img src=\"dot.png\" width=\"15px\" style=\"position:absolute; left:0; top:100px;\"></a>") 

	    mpfile.write("\n")

	    ys[i] = 0.0

	    jfile.write("{")

	    jfile.write("\tChr: ")
	    jfile.write(str(data[i][0]))
	    jfile.write(",\n")

	    jfile.write("\tposition: ")
	    jfile.write(str(data[i][1]))
	    jfile.write(",\n")

	    jfile.write("\tp-value: ")
	    jfile.write(str(data[i][2]))
	    jfile.write(",\n")

	    jfile.write("\tcoordX: ")
	    xxx =  (str(6 + x0p + xs[i] * (x1p - x0p) / MAXxs))

	    jfile.write(xxx)
	    jfile.write(",\n")

	    jfile.write("\tcoordY: ")
	    yyy =  (str( y1p + ys[j] * ((-y1p + y0p) /  (MAXys) ) /1.05  )  ) 

	    jfile.write(y)
#	    jfile.write("\n")

#	    jfile.write("\tSNP:")
#	    jfile.write(str(i))
#	    jfile.write("\n")

#	    jfile.write("\tdata[i][0]: ")	    
#	    jfile.write(data[i][0])
#	    jfile.write(",\n")

#	    jfile.write("\tdata[i][1]: ")
#	    jfile.write(str(data[i][1]))
#	    jfile.write(",\n")

#	    jfile.write("\tdata[i][2]: ")
#	    jfile.write(str(data[i][2]))
#	    jfile.write(",\n")
	    if (j == (int(dotsactive) - 1)):
		jfile.write("\n}")
	    else:
		jfile.write("\n},\n")

#	    mpfile.write("</map>")

    jfile.write("]")
    jfile.close()
    mpfile.close()
#        if ys.aMAXys() > 3:
#    	    print ys
    # plot 0.05 line after multiple testing.

#    ax.axhline(y=-np.log10(0.05 / len(data)), color='0.5', linewidth=2)

    plt.axis('tight')
    plt.xlim(0, xs[-1])
    plt.ylim(ymin=0)
    if ymax is not None: plt.ylim(ymax=ymax)
    plt.xticks([c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], rotation=-90, size=8.5)
    print >>sys.stderr, "saving to: %s" % image_path
    plt.grid()
    plt.savefig(image_path)
    #plt.show()


def get_filehandles(args):
    return (open(a) if a != "-" else sys.stdin for a in args)

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("--no-log", dest="no_log", help="don't do -log10(p) on the value",
            action='store_true', default=False)
    p.add_option("--cols", dest="cols", help="zero-based column indexes to get"
        " chr, position, p-value respectively e.g. %default", default="0,1,2")
    p.add_option("--colors", dest="colors", help="cycle through these colors",
                default="bk")
    p.add_option("--image", dest="image", help="save the image to this file. e.g. %default",
                default="manhattan7.png")
    p.add_option("--title", help="title for the image.", default=None, dest="title")
    p.add_option("--ymax", help="max (logged) y-value for plot", dest="ymax", type='float')
    p.add_option("--sep", help="data separator, default is [tab]",
            default="\t", dest="sep")
    p.add_option("--lines", default=False, dest="lines", action="store_true",
        help="plot the p-values as lines extending from the x-axis rather than"
             " points in space. plotting will take longer with this option.")
    p.add_option("--dotsactive", default="100", dest="dotsactive", help="Number of clickable top p-value")
    p.add_option("--module", default="100", dest="module", help="Module name")
    p.add_option("--path", dest="path", default="./", help="Path where save the results")

    opts, args = p.parse_args()
    if (len(args) == 0):
        sys.exit(not p.print_help())
    fhs = get_filehandles(args)
    columns = map(int, opts.cols.split(","))
#    print ("columns")
#    print (columns)
#    print ("columns ended")
    manhattan(fhs, columns, opts.image, opts.no_log, opts.colors, opts.sep,
            opts.title, opts.lines, opts.ymax, opts.dotsactive, opts.module, opts.path)

if __name__ == "__main__":
    main()