@echo off
chcp 65001 >nul
cd /d "%~dp0"
where julia >nul 2>nul
if errorlevel 1 (
  echo No encuentro julia en el PATH. Instalalo con juliaup o agrega julia al PATH.
  pause
  exit /b 1
)
echo Verificando el codigo...
julia test\runtests.jl
if errorlevel 1 (
  echo Fallaron las pruebas. No se lanza el experimento.
  pause
  exit /b 1
)
echo.
echo ETAPA B: 22 conjeturas abiertas, algoritmo mdes, 10 corridas de 900 s.
echo Se puede cortar con Ctrl+C y volver a lanzar: retoma donde quedo.
echo No dejes que la laptop entre en suspension mientras corre.
echo.
julia -t 8 scripts\experiment.jl --only open --algs mdes --runs 10 --time 900
julia scripts\summary.jl
echo.
echo Listo. Si hubo contraejemplos, estan en results\counterexamples
pause
