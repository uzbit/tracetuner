rem
rem  Build checkbc, checkqv, checkbcphd and checkqvphd on the pc using the GCC compiler
rem

mkdir ..\..\rel
mkdir ..\..\obj
mkdir ..\..\rel\x86-win32
mkdir ..\..\obj\x86-win32
del   ..\..\obj\x86-win32\*.o
copy ..\compute_qv\Btk_atod.h Btk_atod.h
copy ..\compute_qv\Btk_lookup_table.h Btk_lookup_table.h
copy ..\compute_qv\Btk_qv.h Btk_qv.h
copy ..\compute_qv\Btk_qv_data.h Btk_qv_data.h
copy ..\compute_qv\Btk_compute_qv.h Btk_compute_qv.h
copy ..\compute_qv\util.h util.h
copy ..\mktrain\train.h train.h
copy ..\mklut\lut.h lut.h
copy ..\mklut\params.h params.h
copy ..\mktrain\train.h train.h
gcc -D__WIN32 -O3 -c check_data.c -o           ..\..\obj\x86-win32\check_data.o
gcc -D__WIN32 -O3 -c checkbc.c -o              ..\..\obj\x86-win32\checkbc.o
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\checkbc ..\..\obj\x86-win32\check*.o ..\..\lib\x86-win32\*.a
del ..\..\obj\x86-win32\checkbc.o   
gcc -D__WIN32 -O3 -c checkqv.c -o              ..\..\obj\x86-win32\checkqv.o       
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\checkqv ..\..\obj\x86-win32\check*.o ..\..\lib\x86-win32\*.a -lm
del ..\..\obj\x86-win32\checkqv.o
gcc -D__WIN32 -O3 -c checkbcphd.c -o           ..\..\obj\x86-win32\checkbcphd.o
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\checkbcphd ..\..\obj\x86-win32\check*.o ..\..\lib\x86-win32\*.a
del ..\..\obj\x86-win32\checkbcphd.o
gcc -D__WIN32 -O3 -c checkqvphd.c -o           ..\..\obj\x86-win32\checkqvphd.o
gcc -D__WIN32 -O3 -o  ..\..\rel\x86-win32\checkqvphd ..\..\obj\x86-win32\check*.o ..\..\lib\x86-win32\*.a -lm
del ..\..\obj\x86-win32\checkqvphd.o
del Btk_qv.h
del Btk_qv_data.h
del Btk_lookup_table.h
del Btk_compute_qv.h
del util.h
del train.h
del lut.h
del params.h
