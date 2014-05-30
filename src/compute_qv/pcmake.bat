rem
rem  Build ttuner.exe on the pc using the GCC compiler
rem

mkdir ..\..\rel
mkdir ..\..\obj
mkdir ..\..\rel\x86-win32
mkdir ..\..\obj\x86-win32
del  ..\..\obj\x86-win32\*.o
copy ..\mktrain\Btk_compute_match.c Btk_compute_match.c
copy ..\mktrain\Btk_compute_match.h Btk_compute_match.h
copy ..\mktrain\Btk_match_data.h Btk_match_data.h
copy ..\mktrain\Btk_match_data.c Btk_match_data.c
copy ..\mktrain\Btk_sw.c Btk_sw.c
copy ..\mktrain\Btk_sw.h Btk_sw.h
copy ..\mktrain\train.h  train.h
gcc -D__WIN32 -O3 -c ABI_Toolkit.c -o          ..\..\obj\x86-win32\ABI_Toolkit.o
gcc -D__WIN32 -O3 -c Btk_atod.c -o             ..\..\obj\x86-win32\Btk_atod.o
gcc -D__WIN32 -O3 -c Btk_process_raw_data.c -o ..\..\obj\x86-win32\Btk_process_raw_data.o
gcc -D__WIN32 -O3 -c Btk_call_bases.c -o       ..\..\obj\x86-win32\Btk_call_bases.o
gcc -D__WIN32 -O3 -c Btk_compute_qv.c -o       ..\..\obj\x86-win32\Btk_compute_qv.o
gcc -D__WIN32 -O3 -c Btk_compute_tp.c -o       ..\..\obj\x86-win32\Btk_compute_tp.o
gcc -D__WIN32 -O3 -c Btk_compute_tpars.c -o    ..\..\obj\x86-win32\Btk_compute_tpars.o
gcc -D__WIN32 -O3 -c Btk_default_table.c -o    ..\..\obj\x86-win32\Btk_default_table.o
gcc -D__WIN32 -O3 -c Btk_lookup_table.c -o     ..\..\obj\x86-win32\Btk_lookup_table.o
gcc -D__WIN32 -O3 -c Btk_process_peaks.c -o    ..\..\obj\x86-win32\Btk_process_peaks.o
gcc -D__WIN32 -O3 -c Btk_qv_funs.c -o          ..\..\obj\x86-win32\Btk_qv_funs.o
gcc -D__WIN32 -O3 -c Btk_get_mixed_bases.c -o  ..\..\obj\x86-win32\Btk_get_mixed_bases.o
gcc -D__WIN32 -O3 -c Btk_qv_io.c -o            ..\..\obj\x86-win32\Btk_qv_io.o
gcc -D__WIN32 -O3 -c FileHandler.c -o          ..\..\obj\x86-win32\FileHandler.o
gcc -D__WIN32 -O3 -c SCF_Toolkit.c -o          ..\..\obj\x86-win32\SCF_Toolkit.o
gcc -D__WIN32 -O3 -c util.c -o                 ..\..\obj\x86-win32\util.o
gcc -D__WIN32 -O3 -c nr.c   -o                 ..\..\obj\x86-win32\nr.o
gcc -D__WIN32 -O3 -c Btk_match_data.c -o       ..\..\obj\x86-win32\Btk_match_data.o
gcc -D__WIN32 -O3 -c Btk_compute_match.c -o    ..\..\obj\x86-win32\Btk_compute_match.o
gcc -D__WIN32 -O3 -c Btk_sw.c -o               ..\..\obj\x86-win32\Btk_sw.o 
gcc -D__WIN32 -O3 -c context_table.c -o        ..\..\obj\x86-win32\context_table.o
gcc -D__WIN32 -O3 -c tracepoly.c -o            ..\..\obj\x86-win32\tracepoly.o
gcc -D__WIN32 -O3 -c Btk_process_indels.c -o   ..\..\obj\x86-win32\Btk_process_indels.o
gcc -D__WIN32 -O3 -c main.c -o                 ..\..\obj\x86-win32\main.o
gcc -D__WIN32 -O3 -o ..\..\rel\x86-win32\ttuner ..\..\obj\x86-win32\*.o 
del Btk_sw.c 
del Btk_sw.h 
del Btk_match_data.c 
del Btk_match_data.h 
del Btk_compute_match.c 
del Btk_compute_match.h
del ..\..\obj\x86-win32\main.o
