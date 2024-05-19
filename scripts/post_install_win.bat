pip install delvewheel meson ninja

meson setup bbdir
meson compile -C bbdir

REM Repair the wheel using delvewheel
set destDir=%1
set wheel=%2

delvewheel repair --add-path c:/bin;c:/lib;c:/bin/src;c:/lib/src;D:/a/HiGHS/HiGHS/bbdir/src/ -w %destDir% %wheel%
