all: CO.abun abundances.eps abundances.pdf

CO.abun: ../../src/astrochem input.ini
	../../src/astrochem input.ini

abundances.eps: ../../src/plabun CO.abun
	../../src/plabun -o abundances.eps time CO.abun C\(+\).abun C.abun e\(-\).abun

abundances.pdf: abundances.eps
	ps2pdf abundances.eps

CLEANFILES = *.abun *.rout *.eps *.pdf
