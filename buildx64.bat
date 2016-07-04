::This is a batch file for creating a Win64 Visual Studio Project
cmake . -G "Visual Studio 11 Win64"
call "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" x86_amd64
call "C:\Program Files (x86)\Intel\Composer XE 2013\bin\ifortvars.bat" intel64 vs2012
set mypath=%cd%
cd myPath
MSBuild.exe "bftMechanics.sln"