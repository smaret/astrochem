#!@PYTHON@
#
#  plabun - Plot the abundances computed by Astrochem
#
#  Copyright (c) 2006-2011 Sebastien Maret
# 
#  This file is part of Astrochem.
#
#  Astrochem is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  Astrochem is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.

import sys
import getopt
import string
import biggles
from numpy import *

VERSION = "0.2"

def usage():

    # Display usage.

    print """Usage: plabun [options] command file1 [file2...]

Commands:
   time   Plot abundances vs. time in a given shell
   av     Plot abundance vs. av at a given time

Options:
   -h, --help               Display this help
   -V, --version            Display plabun version information
   -o, --output=file        Create a postscript file

   -s, --shell=index        Set the shell number (default 0)
   -t, --time=time          Set the time index (default -1)
   -x, --xrange=xmin,xmax   Set the x axis range
   -y, --yrange=ymin,ymax   Set the y axis range

   -m, --model=file         Specify the model filename
   
See the plabun(1) man page for more information
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

def version():

    # Display version number.

    print "This is plabun, version %s" % VERSION
    print """Copyright (c) 2006-2011 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

def readabun(filename):

    # Read an abund file and return arrays of time and abundances
    
    a = []
    try:
        f = open(filename)
    except:
        sys.stderr.write("plabun: error: can't open %s.\n" % filename)
        sys.exit(1)

    # Skip comments and shell number
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    # Read all lines
    lines = map(string.strip, f.readlines())

    # Construct a list of the elements of the array
    for line in lines:
        newline = []        
        for elem in string.split(line):
	    elem = string.atof(elem)
            newline.append(elem)
        a.append(newline)

    # Create an array from the list
    a = array(a)

    time = a[:,0]
    abund = a[:,1:]

    return time, abund

def readmodel(filename):
    
    # Read a model file (.mdl) and returns the visual extinction and
    # radius (if present in the model file).

    av, radius = [], []
    try:
        f = open(filename)
    except:
        sys.stderr.write("plabun: error: can't open %s.\n" % filename)
        sys.exit(1)

    for line in f.readlines():
        if line[0] == "#":
            continue
        av.append(float(line.split()[1]))
        if len(line.split()) >= 6:
            radius.append(float(line.split()[5]))

    av, radius = array(av), array(radius)
                          
    return av, radius

def speciename(filename):
    
    # Guess the specie name from the filename
    
    specie = "$"
    for char in filename[:-5]:
	if char == '(' or char == ')':
	    continue
        if char == '_':
            break
	elif (char == '+' or char == '-') and char == filename[-7]:
	    specie = specie + "^{" + char + "}"
	elif char.isdigit():
	    specie = specie + "_{" + char + "}"
	else:
	    specie = specie + char
    specie = specie + "$"
    specie = specie.replace("ice", " ice")

    return specie

def main():

    # Parse options and check commands and arguments

    try:
	opts, args = getopt.getopt(sys.argv[1:], "hVo:s:t:x:y:",
				   ["help", "version", "output=",
                                    "shell=", "time=", "xrange=",
                                    "yrange=", "model="])
    except getopt.GetoptError:
	usage()
	sys.exit(1)

    output = None
    s =  0
    t = -1
    xrange = None
    yrange = None
    model = None

    for opt, arg in opts:
	if opt in ("-h", "--help") :
	    usage()
	    sys.exit()
	if opt in ("-V", "--version") :
	    version()
	    sys.exit()
	if opt in ("-o", "--output"):
	    output = arg
	if opt in ("-s", "--shell"):
            try:
                s = int(arg)
            except:
                sys.stderr.write("plabun: error: invalid shell index.\n")
                sys.exit(1)
	if opt in ("-t", "--time"):
            try:
                t = int(arg)
            except:
                sys.stderr.write("plabun: error: invalid time index.\n")
                sys.exit(1)
        if opt in ("-x", "--xrange"):
            try:
                xrange = [float(arg.split(",")[0]), float(arg.split(",")[1])]
            except:
                sys.stderr.write("plabun: error: invalid x axis range.\n")
                sys.exit(1)
        if opt in ("-y", "--yrange"):
            try:
                yrange = [float(arg.split(",")[0]), float(arg.split(",")[1])]
            except:
                sys.stderr.write("plabun: error: invalid y axis range.\n")
                sys.exit(1)
        if opt in ("-m", "--model"):
            model = arg

    if len(args) < 2:
	usage()
	sys.exit(1)
    command = args[0]
    filenames = args[1:]

    if not command in ["time", "av", "radius"]:
	sys.stderr.write("plabun: error: invalid command.\n")
	sys.exit(1)
	
    p = biggles.FramedPlot()
    p.title = "Abundances"
    p.ytitle = "$n(x)/n_{H}$"
    p.ylog = 1
    if xrange:
	p.xrange = xrange
    if yrange:
	p.yrange = yrange
    curves = []

    # Stack for line colors
    # FixMe: use different symbols once we used
    # all colors. Set a maximum number of plots?
    linecolor_stack = ["red", "blue", "green", "yellow", "orange", "cyan"]

    # Read the model file, if given
    if model:
        av, radius = readmodel(model)
    
    # Plot the abundances in each file
    for filename in filenames:
	time, abund = readabun(filename)
	shell = arange(len(time))

	if command == "time":

	    p.xtitle = "Time (yr)"
	    p.xlog = 1

            # Pick-up values that correspond to the given time index
            try:
                abund = abund[:, s]
	    except IndexError:
		sys.stderr.write("plabun: error: shell index is out of bounds.\n")
		sys.exit(1)            

	    # Drop zeros and negative values
            index = abund > 0
	    time = time[index]
	    abund = abund[index]

	    # Check that time and abund contain at least two points
	    if len(time) < 2:
		sys.stderr.write("plabun: warning: %s contains less than two non-zero abundances.\n" % filename)
		continue
        
	    linecolor = linecolor_stack.pop(0)
	    linecolor_stack.append(linecolor)
            c = biggles.Curve(time, abund, linecolor = linecolor, linewidth = 2)
	    c.label = speciename(filename)
	    curves.append(c)
	    p.add(c)

	else:

            # Make sure that the model file was given
            if not(model):
		sys.stderr.write("plabun: error: no model file given.\n")
                sys.exit(1)

            if command == "av":
                p.xtitle = "$A_{v}$ (mag)"
                p.xlog = 0
                xvalues = av
            else:
                # Make sure the radius values are in the model fiel
                if radius.size == 0:
                    sys.stderr.write("plabun: error: no radius in the model file.\n")
                    sys.exit(1)
                p.xtitle = "Radius (AU)"
                p.xlog = 1
                xvalues = radius

            # Pick-up values that correspond to the given time index
            try:
                abund = abund[t, :]
	    except IndexError:
		sys.stderr.write("plabun: error: time index is out of bounds.\n")
		sys.exit(1)            

	    # Drop zeros and negative values
            index = abund > 0
	    xvalues = xvalues[index]
	    abund = abund[index]

	    # Drop zeros and negative values
            index = where(abund > 0)
	    xvalues = xvalues[index]
	    abund = abund[index]

	    # Check that xvalues and abund contain at least two points
	    if len(xvalues) < 2:
		sys.stderr.write("plabun: warning: %s contains less than two non-zero abundances.\n" % filename)
		continue
        
	    linecolor = linecolor_stack.pop(0)
	    linecolor_stack.append(linecolor)
            c = biggles.Curve(xvalues, abund, linecolor = linecolor)
	    c.label = speciename(filename)
	    curves.append(c)
	    p.add(c)
    
    # Draw the plot key
    p.add(biggles.PlotKey(.1, .90, curves))
    if command == "av" or command == "radius":
        p.add(biggles.PlotLabel(.8,.9, "t=%3.1f x 10$^{%i} yr" % 
                                 (time[t]/10**floor(log10(time[t])),
                                  floor(log10(time[t])))))

    # Check that we have at least one abundance to plot
    if len(curves) == 0:
        sys.stderr.write("plabun: error: nothing to plot.\n")
        sys.exit(1)
    
    if output:
	p.write_eps(output)
    else:
        p.show()

main()			  
