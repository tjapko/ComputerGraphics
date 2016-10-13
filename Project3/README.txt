This is the source code and libraries you need in
order to solve the exercises if you are running Windows.

It is recommended but not enforced to use any version
of Microsoft Visual Studio to edit, compile and debug
your code.

The code requires GLUT to compile, which is included
in this package (freeglut-MSVC-2.8.1-1.mp.zip).
The easiest way to make sure Visual Studio knows where
the includes, libaries and binaries are located is by
copying the contents of the /bin, /include and /lib
(in the freeglut directory) to the /bin, /include and /lib
directories found at a path that should be similar to

C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC

If that is not to your liking you will have to link and
include the appropriate files and copy the binaries directly
into your own bin directories.

To compile, run, edit and debug the code, create a new C++ console
application in VS, select the boxes empty project and deselect
precompiled headers  and SDL development lifecycle. Then copy all 
code files into the folder of the solution and in VS right click 
the project and select Add -> Add existing items, and select all .cpp
and .h files within that folder.

The code file that should contain your own code is going
to be myFunctions.cpp, the rest does not need to be touched
but of course you can have a look at everything.
Especially the main.cpp file may help you understand how
to setup an OpenGL viewer, which is surprisingly easy.