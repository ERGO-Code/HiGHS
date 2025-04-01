REM SPDX-FileCopyrightText: 2020 Intel Corporation
REM
REM SPDX-License-Identifier: MIT

set URL=%1

curl.exe --output webimage.exe --url %URL% --retry 5 --retry-delay 5
start /b /wait webimage.exe -s -x -f webimage_extracted --log extract.log
del webimage.exe

webimage_extracted\bootstrapper.exe -s --action install --eula=accept -p=NEED_VS2022_INTEGRATION=0 -p=NEED_VS2019_INTEGRATION=0 --log-dir=.
