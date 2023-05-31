@echo off

set current_dir=%~dp0
echo CPU Architecture is %CPU_ARCHITECTURE% 
echo "copy resource to gefpy"
copy %current_dir%\win\resource\*.dll %current_dir%\gefpy\

python setup.py bdist_wheel

del /a /f /s /q %current_dir%\gefpy\*.dll
del /a /f /s /q %current_dir%\gefpy\*.cpp