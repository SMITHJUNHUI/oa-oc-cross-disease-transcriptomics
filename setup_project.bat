@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup_project.ps1"
if errorlevel 1 (
  echo.
  echo Dependency restore failed.
  pause
  exit /b 1
)
echo.
echo Dependencies restored. You can now run run_project.bat.
pause
