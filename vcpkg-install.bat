REM this file installs the dependencies using VCPKG
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg; .\bootstrap-vcpkg.bat
.\vcpkg.exe integrate install
cd ..
vcpkg\vcpkg install zlib:x64-windows
vcpkg\vcpkg install gsl:x64-windows

setx PATH "%cd%\vcpkg;%PATH%"
set PATH=%cd%\vcpkg;%PATH%