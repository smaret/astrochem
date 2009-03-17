#!@PYTHON@
#
#  plabun - Plot the abundances computed by Astrochem
#
#  Copyright (c) 2006-2009 Sebastien Maret
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

VERSION = "0.1"

# Display usage.

def usage():
    print """Usage: plabun [options] command file1 [file2...]

Commands:
   time   Plot abundances vs. time in a given shell
   av     Plot abundance vs. visual extinction at a given time

Options:
   -h, --help               Display this help
   -V, --version            Display plabun version information
   -o, --output             Create a postscript file

   -s, --shell=index        Set the shell number (default 0)
   -t, --time=time          Set the time
   -x, --xrange=xmin,xmax   Set the x axis range
   -y, --yrange=ymin,ymax   Set the y axis range
   
See the plabun(1) man page for more information
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

# Display version number.

def version():
    print "This is plabun, version %s" % VERSION
    print """Copyright (c) 2006-2009 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

# Read an abund file and return arrays of time, visual extinction 
# and abundances

def readabun(filename):
    
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

    return specie

# Parse options and check commands and arguments
	
def main():

    try:
	opts, args = getopt.getopt(sys.argv[1:], "ho:s:t:x:y:",
				   ["help", "output=", "shell=", 
                                    "time=", "xrange=", "yrange="])
    except getopt.GetoptError:
	usage()
	sys.exit(1)

    output = None
    s =  0
    t = -1
    xrange = None
    yrange = None

    for opt, arg in opts:
	if opt in ("-h", "--help") :
	    usage()
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

    if len(args) < 2:
	usage()
	sys.exit(1)
    command = args[0]
    filenames = args[1:]

    if not command in ["time", "av"]:
	sys.stderr.write("plabun: error: invalid command.\n")
	sys.exit(1)
	
    p = biggles.FramedPlot()
    p.title = "Abundances"
    p.ytitle = "$n(x)/n(H_{2})$"
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
    
    # Plot the abundances in each file
    for filename in filenames:
	time, abund = readabun(filename)
	shell = arange(len(time))

	if command == "time":

	    p.xtitle = "Time"
	    p.xlog = 1

	    # Drop zeros and negative values
	    index = where(abund[:,s] > 0)[0]
	    time = time[index]
	    abund = abund[index]

	    # Check that time and abund contain at least two points
	    if len(time) < 2:
		sys.stderr.write("plabun: warning: %s contains less than two non-zero abundances.\n" % filename)
		continue
        
	    linecolor = linecolor_stack.pop(0)
	    linecolor_stack.append(linecolor)
	    try:
		c = biggles.Curve(time, abund[:, s], linecolor = linecolor, linewidth = 2)
	    except IndexError:
		sys.stderr.write("plabun: error: shell index is out of bounds.\n")
		sys.exit(1)
	    c.label = speciename(filename)
	    curves.append(c)
	    p.add(c)

	else:

	    p.xtitle = "Shell index"
	    p.xlog = 0

	    # Drop zeros and negative values
	    index = where(abund[t,:] > 0)[0]
	    shell = shell[index]
	    abund = abund[index]

	    print abund.shape

	    # Check that time and abund contain at least two points
	    if len(shell) < 2:
		sys.stderr.write("plabun: warning: %s contains less than two non-zero abundances.\n" % filename)
		continue
        
	    linecolor = linecolor_stack.pop(0)
	    linecolor_stack.append(linecolor)
	    try:
		c = biggles.Curve(shell, abund[t,:], linecolor = linecolor)
	    except IndexError:
		sys.stderr.write("plabun: error: time index is out of bounds.\n")
		sys.exit(1)
	    c.label = speciename(filename)
	    curves.append(c)
	    p.add(c)
    
    # Draw the plot key
    p.add(biggles.PlotKey(.1, .90, curves))
    if command == "av":
        p.add(biggles.PlotLabel(.8,.9, "t=%3.1fx10$^{%i} yr" % 
                                 (time[t]/10**log10(time[t]),
                                  floor(log10(time[t])))))

    # Check that we have at least one abundance to plot
    if len(curves) == 0:
        sys.stderr.write("plabun: error: nothing to plot.\n")
        sys.exit(1)
    
    p.show()
    if output:
	p.write_eps(arg)

main()			  
	
	
    
