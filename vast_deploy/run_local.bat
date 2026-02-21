@echo off
REM Run these 5 windows on your local RTX 4060 Ti
REM While Vast.ai handles the other 42

setlocal
set BASE=C:\Users\vasud\nod2-screening-data\fep_pmx\wt_complex

echo ============================================
echo LOCAL FEP RUN - 5 windows on your GPU
echo ============================================
echo.

echo [1/5] Running wt_complex window 00...
cd /d "%BASE%\window_00"
if errorlevel 1 (echo ERROR: Cannot cd to window_00 & goto :error)
conda run -n fep python run_window.py
if errorlevel 1 (echo ERROR: window_00 failed & goto :error)
echo [DONE] wt_complex window 00
echo.

echo [2/5] Running wt_complex window 14...
cd /d "%BASE%\window_14"
if errorlevel 1 (echo ERROR: Cannot cd to window_14 & goto :error)
conda run -n fep python run_window.py
if errorlevel 1 (echo ERROR: window_14 failed & goto :error)
echo [DONE] wt_complex window 14
echo.

echo [3/5] Running wt_complex window 15...
cd /d "%BASE%\window_15"
if errorlevel 1 (echo ERROR: Cannot cd to window_15 & goto :error)
conda run -n fep python run_window.py
if errorlevel 1 (echo ERROR: window_15 failed & goto :error)
echo [DONE] wt_complex window 15
echo.

echo [4/5] Running wt_complex window 16...
cd /d "%BASE%\window_16"
if errorlevel 1 (echo ERROR: Cannot cd to window_16 & goto :error)
conda run -n fep python run_window.py
if errorlevel 1 (echo ERROR: window_16 failed & goto :error)
echo [DONE] wt_complex window 16
echo.

echo [5/5] Running wt_complex window 17...
cd /d "%BASE%\window_17"
if errorlevel 1 (echo ERROR: Cannot cd to window_17 & goto :error)
conda run -n fep python run_window.py
if errorlevel 1 (echo ERROR: window_17 failed & goto :error)
echo [DONE] wt_complex window 17
echo.

echo ============================================
echo === LOCAL GPU COMPLETE ===
echo All 5 windows finished successfully!
echo ============================================
pause
exit /b 0

:error
echo.
echo ============================================
echo === LOCAL RUN FAILED ===
echo Check the error above
echo ============================================
pause
exit /b 1
