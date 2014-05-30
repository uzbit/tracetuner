rem
rem  Build lut.exe on the pc using the GCC compiler
rem

mkdir ..\..\rel
mkdir ..\..\obj
mkdir ..\..\rel\x86-win32
mkdir ..\..\obj\x86-win32
del   ..\..\obj\x86-win32\*.o
copy ..\compute_qv\Btk_qv.h Btk_qv.h
copy ..\compute_qv\util.h util.h
copy ..\compute_qv\Btk_qv_data.h Btk_qv_data.h
copy ..\compute_qv\ABI_Toolkit.h ABI_Toolkit.h
copy ..\compute_qv\FileHandler.h FileHandler.h
copy ..\compute_qv\Btk_get_mixed_bases.h Btk_get_mixed_bases.h
copy ..\compute_qv\Btk_lookup_table.h Btk_lookup_table.h
copy ..\compute_qv\Btk_qv_io.h Btk_qv_io.h
copy ..\compute_qv\context_table.h context_table.h
copy ..\compute_qv\tracepoly.h tracepoly.h
copy ..\compute_qv\Btk_default_table.h Btk_default_table.h
copy ..\compute_qv\Btk_atod.h Btk_atod.h
copy ..\compute_qv\Btk_compute_qv.h Btk_compute_qv.h
copy ..\compute_qv\util.h util.h
copy ..\compute_qv\Btk_qv_funs.h Btk_qv_funs.h
copy ..\compute_qv\nr.h nr.h
copy ..\compute_qv\Btk_process_peaks.h Btk_process_peaks.h
copy ..\compute_qv\Btk_compute_tpars.h Btk_compute_tpars.h
copy ..\compute_qv\Btk_compute_tp.h Btk_compute_tp.h
copy ..\compute_qv\Btk_call_bases.h Btk_call_bases.h
copy ..\compute_qv\Btk_process_raw_data.h Btk_process_raw_data.h
copy ..\mkchk\check_data.h check_data.h
copy ..\mkchk\check_data.c check_data.c
copy ..\mktrain\train.h train.h
gcc -D__WIN32 -O3 -c select.c -o               ..\..\obj\x86-win32\select.o
gcc -D__WIN32 -O3 -c func_name.c -o            ..\..\obj\x86-win32\func_name.o
gcc -D__WIN32 -O3 -c get_thresholds.c -o       ..\..\obj\x86-win32\get_thresholds.o
gcc -D__WIN32 -O3 -c check_data.c -o           ..\..\obj\x86-win32\check_data.o
gcc -D__WIN32 -O3 -c lut.c -o                  ..\..\obj\x86-win32\lut.o
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\lut ..\..\obj\x86-win32\*.o ..\..\lib\x86-win32\*.a -lm
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
del check_data.h
del Btk_qv_funs.h
del Btk_call_bases.h
del Btk_process_raw_data.h
del Btk_process_peaks.h
del Btk_compute_tpars.h
del Btk_compute_tp.h
del nr.h
del train.h
