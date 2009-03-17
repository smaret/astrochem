info_TEXINFOS = astrochem.texi

astrochem_TEXINFOS = intro.texi using.texi examples/input.ini \
	examples/source.mdl examples/CO.abun examples/H3(+).abun \
	examples/HCO(+).abun examples/e.abun examples/abundances.eps \
	examples/abundances.pdf examples/abundances.eps input.texi \
	network.texi examples/network.chm gpl.texi fdl.texi

man_MANS = astrochem.1 plabun.1 plroute.1 chmconvert.1

EXTRA_DIST = $(man_MANS)

$(DESTDIR)$(infodir)/foo.info: foo.info
	$(POST_INSTALL)
         # There may be a newer info file in . than in srcdir.
	-if test -f foo.info; then d=.; \
	          else d=$(srcdir); fi; \
	$(INSTALL_DATA) $$d/foo.info $(DESTDIR)$@; \
	# Run install-info only if it exists.
 	# Use `if' instead of just prepending `-' to the
 	# line so we notice real errors from install-info.
 	# We use `$(SHELL) -c' because some shells do not
 	# fail gracefully when there is an unknown command.
	if $(SHELL) -c 'install-info --version' \
		>/dev/null 2>&1; then \
		install-info --dir-file=$(DESTDIR)$(infodir)/dir \
		$(DESTDIR)$(infodir)/foo.info; \
		else true; fi