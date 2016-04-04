---
project: Q-FEP-MD
project_dir: ../../src
output_dir: ford
project_github: https://github.com/qusers/qsource
project_website: http://github.com 
summary: Developers documentation for Q
author: Johan Åqvist and others.
author_description: Molecular Simulation Lab. Chief
github: https://github.com/esguerra
email: mauricio.esguerra@gmail.com
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         protected
         private
source: true
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: 
json_module: http://jacobwilliams.github.io/json-fortran/ 
futility: http://cmacmackin.github.io 
license: by-nc 
extra_filetypes: sh #
md_extensions: markdown.extensions.toc(anchorlink=True)
               markdown.extensions.smarty(smart_quotes=False)
---
--------------------

[TOC]

Brief description
-----------------

A user-friendly and object-oriented API for reading and writing JSON files, written in
modern Fortran (Fortran 2003+).  The source code is
[a single Fortran module](|url|/module/json_module.html) file
([json_module.F90](|url|/sourcefile/json_module.f90.html)).

License
-------

The JSON-Fortran source code and related files and documentation are
distributed under a permissive free software license (BSD-style).  See
the
[LICENSE](http://jacobwilliams.github.io/json-fortran/page/development-resources/LICENSE.html)
file for more details.
