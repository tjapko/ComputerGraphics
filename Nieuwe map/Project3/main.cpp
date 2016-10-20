/** 
This file contains everything to set up an OpenGL mesh viewer,
but for your exercises you do not need to touch it.
(But you can, if you want to try out stuff).
myFunctions.cpp contains the blank functions that
you should implement.
*/

#include <GL/glut.h>
#include "traqueboule.h"
#include "mesh.h"
#include "myFunctions.h"

static Mesh myMesh;
bool wireFrameMode = true;

/* Function that gets called
 from the OpenGL main loop
 everytime the screen has to be redrawn. */
void display(void) {
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
    tbVisuTransform();
    draw( myMesh );    
    glutSwapBuffers();
}

/* Function that gets called
 from the OpenGL main loop
once the window is being resized. */
void reshape(int w, int h) {
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective (50, (float)w/h, 1, 10);
    glMatrixMode(GL_MODELVIEW);
}

/* Function that gets called
 from the OpenGL main loop
 once a keypress is being registered.
 This connects the keys 1-5 to the
 functions in myFunctions.cpp*/
void keyboard(unsigned char key, int x, int y) {		
	switch (key) {
	case 'w':
		if (wireFrameMode) {
			wireFrameMode = false;
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		else {
			wireFrameMode = true;
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		glutPostRedisplay();
		break;
	case '1':
		myFunction1(myMesh);
		break;
	case '2':
		myFunction2(myMesh);
		break;
	case '3':
		myFunction3(myMesh);
		break;
	case '4':
		myFunction4(myMesh);
		break;
	case '5':
		myFunction5(myMesh);
		break;
	case '6':
		myFunction6(myMesh);
		break;
	case '7':
		myFunction7(myMesh);
		break;
	case '8':
		myFunction8(myMesh);
		break;
	case '9':
		myMesh = myFunction9(myMesh);
		break;
	case '0':
		myFunction0(myMesh);
		break;
	case 's':
		myMesh = myFunctionS(myMesh);
		break;
	case 'h':
		myFunctionH();
		break;
	}
}

/** Initializes an OpenGL window, loads the mesh,
 sets up the camera and mouse navigation functionality. */
void init(int argc, char ** argv) {
	glutInit(&argc, argv);

	myMesh.loadMesh(argc == 2 ? argv[1] : "spikeyBall.obj");

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	glutInitWindowPosition(600, 100);
    glutInitWindowSize(800,800);
	glutCreateWindow("Model Viewer");

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0,0,-4);
    tbInitTransform();
    tbHelp();

    glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMouseFunc(tbMouseFunc);
	glutMotionFunc(tbMotionFunc);
	glutKeyboardFunc(keyboard);

    glEnable( GL_DEPTH_TEST );
    glShadeModel(GL_SMOOTH);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT);

	if (wireFrameMode) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	/*glEnable(GL_LIGHTING);	glEnable(GL_LIGHT0);	GLfloat lightPos[] = { 0, 1.0, 1.0, 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);*/
}

/** Initialize and then go straight
 to the OpenGL main loop. 
 From there the functions display,
 keyboard, reshape, ... are called.*/
int main(int argc, char **argv) {

	init(argc, argv);

	glutMainLoop();
	return 0;
}