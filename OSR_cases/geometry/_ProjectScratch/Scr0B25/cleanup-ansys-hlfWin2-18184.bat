@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 17964)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 14540)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 15088)
if /i "%LOCALHOST%"=="hlfWin2" (taskkill /f /pid 18184)

del /F cleanup-ansys-hlfWin2-18184.bat
