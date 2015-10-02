def citation(bibtex=False):
    """Return citation string. Use bibtex=True to get
    a bibtex-formatted entry.

    Print the return value to see line breaks."""

    s = """Fidimag - FInite DIfference microMAGnetic simulator. 
https://github.com/fangohr/fidimag
(C) University of Southampton, 2012 - 2015"""

    s_bib = """@Misc{Fidimag,
  author = 	 {University of Southampton},
  title = 	 {Fidimag - FInite DIfference microMAGnetic simulator},
  year = 	 {2015},
  note = 	 {https://github.com/fangohr/fidimag}
}"""
    
    if bibtex: 
        return s_bib
    else:
        return s


