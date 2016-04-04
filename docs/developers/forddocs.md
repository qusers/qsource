---
project: Q Developers
project_dir: ../../src
output_dir: doc
project_github: https://github.com/jacobwilliams/json-fortran
project_download: https://github.com/jacobwilliams/json-fortran/releases/latest
summary: Q. A Free Energy Simulation Machine
author: Johan Aqvist
github: https://github.com/jacobwilliams
website: http://qdyn.no-ip.org
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         protected
         private
source: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
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