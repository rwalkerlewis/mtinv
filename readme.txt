MTINV Version 4.0.1
Sun Feb  4 21:08:15 PST 2024

Updated(see bottom)
Tue Sep 19 00:10:53 PDT 2023
Sun Feb  4 21:08:22 PST 2024

Dependencies:
Requires: GCC and GFortran compilers
Optional: sqlite3, GMT version +4.5.x, GMT +5.x.x, GMT -6.3.x

Tested Systems: 
	1. MacOS Ventura intel x86_64  gcc 11.2   20211101
	2. MacOS Ventura Apple-m2      gcc 13.0.1 20230202 (experimental)
	3. Redhat Linux intel x86_64   gcc 4.8.5  20150623 (Red Hat 4.8.5-44)
	4. Linux Mint                  gcc 11.4.0

Instructions: Type make clean; make all

Executables
FlinnEngdahl:                Mach-O 64-bit executable arm64
area_sphere:                 Mach-O 64-bit executable arm64
clean_exec_in_directory.csh: C shell script text executable, ASCII text, with very long lines (316)
glib2inv:                    Mach-O 64-bit executable arm64
glib2inv_join:               Mach-O 64-bit executable arm64
grn2Mxy:                     Mach-O 64-bit executable arm64
grn2db:                      Mach-O 64-bit executable arm64
grnlib2sac:                  Mach-O 64-bit executable arm64
list_MTdb.csh:               C shell script text executable, ASCII text
makepar:                     Mach-O 64-bit executable arm64
mkgrnlib:                    Mach-O 64-bit executable arm64
mtbestfit:                   Mach-O 64-bit executable arm64
mteig:                       Mach-O 64-bit executable arm64
mtinv:                       Mach-O 64-bit executable arm64
mtscreen.py:                 Python script text executable, Unicode text, UTF-8 text
multithread_mkgrnlib:        Mach-O 64-bit executable arm64
pltmod:                      Mach-O 64-bit executable arm64
print_MTdb.csh:              C shell script text executable, ASCII text
remove_MTdb.csh:             C shell script text executable, ASCII text
rename_SACPZs.csh:           C shell script text executable, ASCII text
renamesac:                   Mach-O 64-bit executable arm64
sac2gmtmap:                  Mach-O 64-bit executable arm64
sac2xy:                      Mach-O 64-bit executable arm64
sacdata2inv:                 Mach-O 64-bit executable arm64
sacmerge:                    Mach-O 64-bit executable arm64
sacqc:                       Mach-O 64-bit executable arm64
scripts_original:            directory
setupMT:                     Mach-O 64-bit executable arm64
stats:                       Mach-O 64-bit executable arm64
switch_color:                Mach-O 64-bit executable arm64
unpack.csh:                  C shell script text executable, ASCII text
updateMTdb:                  Mach-O 64-bit executable arm64
whatshere:                   Mach-O 64-bit executable arm64


Do not forget to set ./mtinv.4.x.x/bin directory to your executable PATH 
shell variables.  See file environmental_variables.csh for additional 
required and optional environmental variables.

Bugfixes
Tue Sep 19 00:12:47 PDT 2023
1. sacdata2inv - memory leak fix with linux mint gcc compilers
2. sacdata2inv - renamed fmul() to scale_data(), wierd problem where linux mint gcc compiler failed to reconize the subroutine
3. added GMT version 4.5.x support to mteig plotting
4. upgraded  GMT version 4.5.x support to mtinv plotting waveform plot
5. various sprintf and strcpy issues when compiling with linux mint gcc compilers

New features:
Sun Feb  4 21:08:44 PST 2024
1. updates to setupMT -hspeec96
2. hspec96_to_grnlib.c


-- Gene Ichinose
