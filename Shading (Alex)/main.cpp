#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "traqueboule.h"
#include "mesh.h"
#include "yourcode.h"

using namespace std;
void computeLighting();
void dealWithUserInput(int x, int y);
unsigned int W_fen = 800;  // largeur fenetre
unsigned int H_fen = 800;  // hauteur fenetre


//________________________________
//________________________________
//________________________________
//________________________________
//________________________________
//Start reading here!!!

//Background color
float BackgroundColor[]={0.2,0.2,0};
// Different display modes
bool Debug=false;
bool DiffuseLighting=false;
bool PhongSpecularLighting=false;
bool BlinnPhongSpecularLighting=false;
bool ToonLightingDiffuse=false;
bool ToonLightingSpecular=false;

enum InterfaceLightPlacementValues{LIGHT_SPHERE_PLACEMENT=0, LIGHT_SHADOW_PLACEMENT, LIGHT_SPECULAR_PLACEMENT};
InterfaceLightPlacementValues InterfaceLightPlacement=LIGHT_SPHERE_PLACEMENT;

//There are several light positions and you can later add light sources
//for the moment, just work with LightPos[0]
//to set the light position to the camera location, you can use 'l'
std::vector<Vec3Df> LightPos;
//Later in the exercise, you can also modify and use a special color per light source.
//The LightColor array is used for this purpose:
std::vector<Vec3Df> LightColor;
//the current light
int SelectedLight=0;

//The position of the camera
//this variable will be updated by the program! Do not change it manually!
//Of course, you can USE its value
Vec3Df CamPos = Vec3Df(0.0f,0.0f,-4.0f);


//Pressing 's' will display the currently chosen vertex
bool ShowSelectedVertex=false;

//The index of the chosen vertex is stored in this variable
//the value is -1 if no vertex is chosen
//to choose a vertex, hover over it and press space
int SelectedVertex=-1;

//per vertex attributes, useful for materials - see later exercises
std::vector<Vec3Df> MeshMaterial;

//this is the MAIN function in which you have to work
//it is automatically iterated over the entire mesh and all light sources (light is the current light index)
//it receives the current vertex location, its normal, and the index of the current vertex (yes, this is somewhat redundant because the mesh  
Vec3Df computeLighting(Vec3Df & vertexPos, Vec3Df & normal, unsigned int light, unsigned int index)
{
	Vec3Df result(0,0,0);
	//do not change any global variables here! This function will be called for EACH vertex 
	//of the mesh, so your change would happen several times
	if (Debug)
	{
		return debugColor(index);
	}
	
	if (DiffuseLighting)
		{
			result+=diffuseOnly(vertexPos, normal, LightPos[light], index);
		}
	if (PhongSpecularLighting)
		{
			result+=phongSpecularOnly(vertexPos, normal, LightPos[light], CamPos, index);
		}
	else if (BlinnPhongSpecularLighting)
		{
			result+=blinnPhongSpecularOnly(vertexPos, normal, LightPos[light], CamPos, index);
		}
	if (ToonLightingDiffuse)
		{
			//overwrite previous stuff
			result=toonShadingNoSpecular(vertexPos, normal, LightPos[light], index);
		}
	if (ToonLightingSpecular&&!DiffuseLighting&&!PhongSpecularLighting&&!BlinnPhongSpecularLighting)
	{
		result+=toonShadingOnlySpecular(vertexPos, normal, LightPos[light], CamPos, index);
	}
	return result;
}

//User interaction - when the user chooses a vertex, you receive its position, normal, its index 
//you can use it to NOW modify all global variables, such as the light position, or change material properties etc.
void userInteraction(const Vec3Df & selectedPos, const Vec3Df & selectedNormal, int selectedIndex)
{
	switch(InterfaceLightPlacement)
	{
	case LIGHT_SPHERE_PLACEMENT:
		{
			LightPos[SelectedLight]=userInteractionSphere(selectedPos, CamPos);
			break;
		}
	case LIGHT_SHADOW_PLACEMENT:
		{
			LightPos[SelectedLight]=userInteractionShadow(selectedPos, selectedNormal,LightPos[SelectedLight]);
			float condition = Vec3Df::dotProduct(LightPos[SelectedLight], selectedNormal);
			cout << condition;
			break;
		}
	case LIGHT_SPECULAR_PLACEMENT:
		{
			LightPos[SelectedLight]=userInteractionSpecular(selectedPos, selectedNormal,LightPos[SelectedLight], CamPos);
			break;
		}
	}
}


// prise en compte du clavier
//Vous allez ajouter quelques fonctionalites pendant le TP
//ce qui est important pour vous: key contient le caractère correspondant à la touche appuyé par l'utilisateur

void keyboard(unsigned char key, int x, int y)
{		
	switch (key)
    {
	case '0':
			Debug=!Debug;
			break;

	case '1': 
		DiffuseLighting=!DiffuseLighting;
		break;

	case '2': 
		PhongSpecularLighting=!PhongSpecularLighting;
		break;

	case '3': 
		BlinnPhongSpecularLighting=!BlinnPhongSpecularLighting;
		break;

	case '4': 
		ToonLightingDiffuse=!ToonLightingDiffuse;
		break;

	case '5': 
		ToonLightingSpecular=!ToonLightingSpecular;
		if (ToonLightingSpecular)
			ToonLightingDiffuse=true;
		break;

	case '6':
		cout<<"Number keys from 6 on not used"<<endl<<endl;
		return;
	case 'M':
	case 'm':
		{
			int temp=InterfaceLightPlacement;
			++temp;
			if (temp>2)
				temp=0;
			InterfaceLightPlacement=(InterfaceLightPlacementValues) temp;
			break;
		}
	case 'l':
		{
			LightPos[SelectedLight]=getCameraPosition();
			return;
		}
	case 'L':
		{
			LightPos.push_back(getCameraPosition());
			LightColor.push_back(Vec3Df(1,1,1));
			return;
		}
	case '+':
		{
			++SelectedLight;
			if (SelectedLight>=(int)LightPos.size())
				SelectedLight=0;
			return;
		}
	case '-':
		{
			--SelectedLight;
			if (SelectedLight<0)
				SelectedLight=LightPos.size()-1;
			return;
		}
	case 'N':
		{	//reset all lights
			LightPos.resize(1);
			LightPos[0]=Vec3Df(0,0,3);
			LightColor.resize(1);
			LightColor[0]=Vec3Df(1,1,1);
			SelectedLight=0;
			return;
		}

	case 's':
		{
			ShowSelectedVertex=!ShowSelectedVertex;
			return;
		}
	case ' ':
		{
			//You do not need to look at the function below
			//it does some computations then calls dealWithUserInput
			dealWithUserInput(x,y);
			return;
		}
	default:
		yourKeyboardFunction(key);
	}

	if (Debug)
	{
		cout<<endl<<("DEBUG MODE - only direct color function is called")<<endl;
		return;
	}
	if (!ToonLightingDiffuse)
	{
		cout<<"REALISTIC SHADING!"<<endl;
		cout<<"__________________"<<endl;
		if (DiffuseLighting)
		{
			cout<<("DiffuseLighting ON")<<endl;
		}
		else
		{
			cout<<("DiffuseLighting OFF")<<endl;
		}

		if (PhongSpecularLighting)
		{
			cout<<"PhongSpecularLighting ON"<<endl;
			cout<<"BlinnPhongSpecularLighting IGNORED"<<endl;
		}
		else if (BlinnPhongSpecularLighting)
		{
			cout<<"PhongSpecularLighting OFF"<<endl;
			cout<<"BlinnPhongSpecularLighting ON"<<endl;
		}
		else
		{
			cout<<"PhongSpecularLighting OFF"<<endl;
			cout<<"BlinnPhongSpecularLighting OFF"<<endl;
		}
	}
	else 
	{
		cout<<"TOON SHADING!"<<endl;
		cout<<"_____________"<<endl;
		if (ToonLightingDiffuse)
		{
			cout<<"ToonLightingDiffuse ON"<<endl;
		}
		else
		{
			cout<<"ToonLightingDiffuse OFF"<<endl;
		}

		if (ToonLightingSpecular)
		{
			cout<<"ToonLightingSpecular ON"<<endl;
		}
		else
		{
			cout<<"ToonLightingSpecular OFF"<<endl;
		}
	}
	if (InterfaceLightPlacement==LIGHT_SPHERE_PLACEMENT)
	{
		cout<<"Interaction: LIGHT_SPHERE_PLACEMENT"<<endl;
	}
	if (InterfaceLightPlacement==LIGHT_SHADOW_PLACEMENT)
	{
		cout<<"Interaction: LIGHT_SHADOW_PLACEMENT"<<endl;
	}
	if (InterfaceLightPlacement==LIGHT_SPECULAR_PLACEMENT)
	{
		cout<<"Interaction: LIGHT_SPECULAR_PLACEMENT"<<endl;
	}

}

//Everything below, you do NOT have to read!!!
//________________________________
//________________________________
//________________________________
//________________________________
//________________________________


//Lighting of the model (you should not need to touch this one...
std::vector<Vec3Df> lighting;

/************************************************************
 * Fonction pour initialiser le maillage
 ************************************************************/
void init(const char * fileName){

	//this function loads a mesh
    MyMesh.loadMesh(fileName);
	lighting.resize(MyMesh.vertices.size());
	initStudentVariables();

	MeshMaterial.resize(MyMesh.vertices.size());
	for (unsigned int i=0; i<MyMesh.vertices.size();++i)
		MeshMaterial[i]=Vec3Df(0,0,0);
		
	LightPos.push_back(Vec3Df(0,0,3));
	LightColor.push_back(Vec3Df(1,1,1));
	computeLighting();
}



/************************************************************
 * Appel des différentes fonctions de dessin
************************************************************/


void dealWithUserInput(int x, int y)
{
	Vec3Df worldPoint=getWorldPositionOfPixel(x, H_fen-y);
	SelectedVertex=MyMesh.getClosestVertexIndex(CamPos, worldPoint-CamPos);
	if (SelectedVertex>=0)
	{
		Vec3Df selectedPos=MyMesh.vertices[SelectedVertex].p;
		Vec3Df selectedNormal=MyMesh.vertices[SelectedVertex].n;
		userInteraction(selectedPos, selectedNormal, SelectedVertex);				
	}
}

void dessiner( )
{

	glPointSize(10);
	glBegin(GL_POINTS);
	//LightColor
	glColor3f(1,0,0);	
	
	for (unsigned int i=0; i<LightPos.size();++i)	
	{	
		glVertex3f(LightPos[i][0],LightPos[i][1],LightPos[i][2]);
	}
	glEnd();
	
	glPointSize(40);
	glColor3f(1,1,0);	
	glBegin(GL_POINTS);
	glVertex3f(LightPos[SelectedLight][0],LightPos[SelectedLight][1],LightPos[SelectedLight][2]);
	glEnd();

	Vec3Df p;
	if (ShowSelectedVertex&&SelectedVertex>=0)
	{
		p=MyMesh.vertices[SelectedVertex].p;
		glBegin(GL_POINTS);
		glVertex3f(p[0],p[1],p[2]);
		glEnd();
	}

	MyMesh.drawWithColors(lighting);
}




void idle()
{
	CamPos=getCameraPosition();

	computeLighting();

	glutPostRedisplay();
}

void display(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);


void computeLighting()
{
	std::vector<Vec3Df> *result=&lighting;


	for (unsigned int i=0; i<MyMesh.vertices.size();++i)
	{
		(*result)[i]=Vec3Df();
		for (unsigned int l=0; l<LightPos.size();++l)
			(*result)[i]+=computeLighting(MyMesh.vertices[i].p, MyMesh.vertices[i].n, l, i);
	}
}



/************************************************************
 * Programme principal
 ************************************************************/
int main(int argc, char** argv)
{
    glutInit (&argc, argv);

    init(argc == 2 ? argv[1] : "bunny.obj");

    // couches du framebuffer utilisees par l'application
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

    // position et taille de la fenetre
    glutInitWindowPosition(200, 100);
    glutInitWindowSize(W_fen,H_fen);
    glutCreateWindow(argv[0]);	

    // Initialisation du point de vue
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0,0,-4);
    tbInitTransform();     // initialisation du point de vue
    tbHelp();                      // affiche l'aide sur la traqueboule

    glDisable( GL_LIGHTING );
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    
    // cablage des callback
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutMouseFunc(tbMouseFunc);    // traqueboule utilise la souris
    glutMotionFunc(tbMotionFunc);  // traqueboule utilise la souris
    glutIdleFunc(idle);


    // Details sur le mode de tracé
    glEnable( GL_DEPTH_TEST );            // effectuer le test de profondeur
    glShadeModel(GL_SMOOTH);

    // Effacer tout
    glClearColor (BackgroundColor[0],BackgroundColor[1], BackgroundColor[2], 0.0);
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT); // la couleur et le z
   

	cout<<"Program Usage:"<<endl;
	cout<<"0 - activate Debug"<<endl;
	cout<<"__________________"<<endl;
	cout<<"1 - Diffuse Lighting on"<<endl;
	cout<<"2 - Phong Specularities"<<endl;
	cout<<"3 - Blinn-Phong Specularities"<<endl;
	cout<<"4 - Toon-Shading"<<endl;
	cout<<"5 - Toon Specularities"<<endl;
	cout<<"______________________"<<endl;
	cout<<"m - Change Interaction Mode - to influence light source"<<endl;
	cout<<"l - place the light source at the current camera position"<<endl;
	cout<<"L - add an additional light source"<<endl;
	cout<<"+ - choose next light source"<<endl;
	cout<<"- - choose previous light source"<<endl;
	cout<<"N - clear all light sources and reinitialize with one"<<endl;
	cout<<"s - show selected vertices"<<endl;
	cout<<"SPACE - replaces mouse click for selection, will then call your light placement function"<<endl;

    // lancement de la boucle principale
    glutMainLoop();

    return 0;  // instruction jamais exécutée
}


/************************************************************
 * Fonctions de gestion opengl à ne pas toucher
 ************************************************************/
// Actions d'affichage
// Ne pas changer
void display(void)
{
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT); // la couleur et le z

	glLoadIdentity();  // repere camera

    tbVisuTransform(); // origine et orientation de la scene

    dessiner( );    

    glutSwapBuffers();
}

// pour changement de taille ou desiconification
void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective (50, (float)w/h, 1, 10);
    glMatrixMode(GL_MODELVIEW);
}

