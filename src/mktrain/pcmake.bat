rem
rem  Build train.exe and trainphd.exe on the pc using the GCC compiler
rem

mkdir ..\..\rel
mkdir ..\..\obj
mkdir ..\..\rel\x86-win32
mkdir ..\..\obj\x86-win32
copy ..\compute_qv\Btk_qv.h Btk_qv.h
copy ..\compute_qv\util.h util.h
copy ..\compute_qv\Btk_qv_data.h Btk_qv_data.h
copy ..\compute_qv\ABI_Toolkit.h ABI_Toolkit.h
copy ..\compute_qv\FileHandler.h FileHandler.h
copy ..\compute_qv\Btk_compute_tpars.h Btk_compute_tpars.h
copy ..\compute_qv\Btk_get_mixed_bases.h Btk_get_mixed_bases.h
copy ..\compute_qv\Btk_lookup_table.h Btk_lookup_table.h
copy ..\compute_qv\Btk_qv_io.h Btk_qv_io.h
copy ..\compute_qv\context_table.h context_table.h
copy ..\compute_qv\tracepoly.h tracepoly.h
copy ..\compute_qv\Btk_default_table.h Btk_default_table.h
copy ..\compute_qv\Btk_atod.h Btk_atod.h
copy ..\compute_qv\Btk_compute_qv.h Btk_compute_qv.h
gcc -D__WIN32 -O3 -c train_data.c -o           ..\..\obj\x86-win32\train_data.o
gcc -D__WIN32 -O3 -c Btk_compute_match.c -o    ..\..\obj\x86-win32\Btk_compute_match.o
gcc -D__WIN32 -O3 -c Btk_match_data.c -o       ..\..\obj\x86-win32\Btk_match_data.o
gcc -D__WIN32 -O3 -c Btk_sw.c -o               ..\..\obj\x86-win32\Btk_sw.o
gcc -D__WIN32 -O3 -c train.c -o                ..\..\obj\x86-win32\train.o
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\train ..\..\obj\x86-win32\*.o 
del ..\..\obj\x86-win32\train.o
gcc -D__WIN32 -O3 -c trainphd.c -o             ..\..\obj\x86-win32\trainphd.o
gcc -D__WIN32 -O3 -o ..\..\rel\x86-win32\trainphd ..\..\obj\x86-win32\*.o 
del ..\..\obj\x86-win32\trainphd.o
del Btk_qv.h
del util.h
del Btk_qv_data.h
del ABI_Toolkit.h
del FileHandler.h
del Btk_compute_tpars.h
del Btk_get_mixed_bases.h
del Btk_lookup_table.h
del Btk_qv_io.h
del tracepoly.h
del context_table.h
