#Q Maintenance guide

####John Marelius, August 29, 2000*

#####Update:Mauricio Esguerra March 13, 2014*


This document describes how to maintain the Q programs at the Åqvist
group.  Recently the code, documentation, and some scripts have been
added to version control using github, these documents still need some
updating. 


##File locations

The current source code for Q and documentation has been pushed to a
github private organization located at https://github.com/qusers


+--------------------------------------+--------------------------------------+
| **What**                             | **Where**                            |
+--------------------------------------+--------------------------------------+
| source code, makefile, history files | qsource/src /  history/              |
+--------------------------------------+--------------------------------------+
| manual                               | qsource/documentation/manual/qman5.p |
|                                      | df                                   |
+--------------------------------------+--------------------------------------+
| license agreement                    | qsource/documentation/license.pdf    |
+--------------------------------------+--------------------------------------+
| this document                        | qsource/documentation/qmaintenancegu |
|                                      | ide.docx                             |
+--------------------------------------+--------------------------------------+
| force field files                    | qsource/ff/                          |
+--------------------------------------+--------------------------------------+
| scripts                              | scripts/                             |
+--------------------------------------+--------------------------------------+
| tests                                | qsource/tests/                       |
+--------------------------------------+--------------------------------------+
| Q web site                           | pending                              |
+--------------------------------------+--------------------------------------+



##Updating Q

You need to be a member of the owners team at the github repository. 

Updating source code on the web server

NOTE CHANGE THIS TO GIT WAY

After modifying any of the master source files in g:\\src\\q4:

-   Update the version number and last-modified-date of the file in
    question, e.g. md.f90, and in the main program, e.g. qdyn.f90.
-   Run the script g:\\src\\tarQ4.csh to update the files on the web
    server.

Updating force field files

NOTE CHANGE THIS TO GIT WAY

After modifying the master files in g:\\FF4, run the script
g:\\FF4\\ff2web.csh to update the files on the web server (both
separated files and compressed archive files).

Updating the manual

NOTE CHANGE THIS TO GIT WAY

After making changes to the manual, update the PDF version on the web
server as follows:

-   Log on to Gem (outside Erling’s office).
-   Open the manual (g:\\doc\\Qman4.doc)
-   Print it to the print driver called ‘Acrobat Distiller’. Check the
    options to see that Letter size paper is used (so that our US users
    can print without problem).
-   When the Adobe Acrobat window appears, save the PDF file as
    Qman4.pdf in the wwwroot\\Q directory on the web server.

##Building executables

This section describes how to transfer source code, compile and package
the executables.

###Windows

-   Not supported for now, we need a developer here.

\

###Mac OSX (mavericks)

\

###Linux CentOS 6.5 (triolith)

\

###Linux (csb)

-   Transfer the whole source code archive (update it first!) by ftp to
    a directory named q4 in your home directory on genome.ibg.uu.se.
-   Log on to genome.ibg.uu.se with SSH.
-   Unpack the source code archive:\
     gunzip -c Q4\_src.tar.gz | tar -xvf -

-   Alternatively, replace steps 1 and 2 by transferring individual
    source files.

-   Remove old object files:\
     make clean
-   Build the executables for use on MIPS R10000 processors:\
     make irix6.4
-   Package the executables:\
     ./package\_R10k.csh (copy this script from ˜john/src/q4)
-   Remove object files\
     make clean
-   Build the executables for use on MIPS R4000 processors:\
     make irix6.2
-   Package the executables:\
     ./package\_R4k.csh (copy this script from ˜john/src/q4)
-   Retrieve the archive files Q4\_IRIX6.4\_R10k.tar.gz and
    Q4\_IRIX6.2\_R4k.tar.gz by ftp to the directory
    wwwroot\\Q\\download\\bin on the web server.

###Linux (abisko)

\

###Linux (tintin)

\

##Updating the Q web site

The Åqvist group web site http://aqvist.bmc.uu.se and the Q web site
http://aqvist.bmc.uu.se/Q can both be edited using VisualPage. The
VisualPage project file is w:\\aqvist\_dev.vpp. The HTML and other files
in this project are all in the w:\\wwwdev directory. Edit the files
there and then copy the modified files, or all files, to w:\\wwwroot
which is the actual web server directory.

The E-mail list

The home page of the list is https://groups.google.com/d/forum/qmoldyn.
To see the member list etc, log in to your google account, you should be
a group owner, go to the group web page and at the upper right side of
the page click on Manage. This will show you all the users at the group
(which is configured as a mailing list), and also a lot of configuration
information for the list.

\

To send a message to the list, address it to qmoldyn@googlegroups.com