#include <vector>
#include "mesh.h"

//THIS IS THE ONLY FILE YOU NEED TO MODIFY and SUBMIT!!! NOTHING ELSE!!!
//Please RENAME as described below and then send it to assignments.eisemann@gmail.com 
//(I will not use this email for any other purposes than collecting your files)
//You should send per student only ONE file (this one).
//Each file should be renamed using your studentids.
//E.g., Imagine two people worked together and their ids are 12 and 14, 
//then you should send two files called "12.h" and "14.h"
//Any deviation from this naming convention will lead to 0 points because the exercise is verified automatically!!!
//Good Luck!

using namespace std;
//Global array to store mesh material properties and algorithmic parameters
vector<Vec3Df> Kd;//diffuse coefficient per vertex
vector<Vec3Df> Ks;//specularity coefficient per vertex
vector<float> Shininess;//exponent for phong and blinn-phong specularities
int ToonDiscretize=4;//number of levels in toon shading
float ToonSpecularThreshold=0.49;//threshold for specularity

//Mesh - will be filled and loaded outside.
Mesh MyMesh;


//Helper function that you can ignore!
void initStudentVariables()
{
	//this is needed so that your table has the same size as the number of vertices.
	//later, the vertex indices received in the functions correspond to the same location in your vector.
	//in other words, you can store some information per vertex here.
	Kd.resize(MyMesh.vertices.size(), Vec3Df(0.5,0.5,0.5));
	Ks.resize(MyMesh.vertices.size(), Vec3Df(0.5,0.5,0.5));
	Shininess.resize(MyMesh.vertices.size(), 3);
}


//for debugging purposes or variable changes (e.g., modify the toon threshold as done below)
//please notice that some keys are already in use!
void yourKeyboardFunction(unsigned char key)
{
	cout<<"Key not used so far, so executing your code!"<<endl;
	
	//recolor the mesh 
	switch(key)
	{
		case 't': 
			ToonSpecularThreshold-=0.001;
		break;
		case 'T': 
			ToonSpecularThreshold+=0.001;
		break;
		case 'd': 
			ToonDiscretize-=1;
			if (ToonDiscretize<=0)
				ToonDiscretize=1;
		break;
		case 'D': 
			ToonDiscretize+=1;
		break;
		
		//in case you want to add colors! - Not mandatory

		case 'r': //decrase diffuse Kd coefficient in the red channel by 0.01
		break;
		case 'R': //increase diffuse Kd coefficient in the red channel by 0.01
		break;
		case 'g'://same for green
		break;
		case 'G':
		break;
		case 'b'://same for blue
		break;
		case 'B':
		break;
	}
	
	cout<<"ToonSpecular"<<ToonSpecularThreshold<<endl;
	cout<<"Toon Discretization levels"<<ToonDiscretize<<endl;

}


//Debug function
Vec3Df debugColor(unsigned int index)
{	//this function you can use in any way you like!
	//e.g., for debugging purposes!
	return MyMesh.vertices[index].n;
	//or 
	//return Kd[index];
}


///////////////
///Shading
///////////////
//standard lambertian shading: Kd * dot(N,L), clamped to zero when negative. Where L is the light vector
//
Vec3Df diffuseOnly(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, unsigned int index)
{	
	Vec3Df Kd = dot
	Vec3Df D = Id * Kd * cos(theta);

	return Vec3Df(1,0,0);
}


//Phong (!) Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
//Follow the course, only calculate Ks pow(dot(V,R),shininess), where V is the view vector and R is the Reflection vector of the light (like in pool billard computed from the LightPos, vertexPos and normal).
//When computing specularities like this, verify that the light is on the right side of the surface, with respect to the normal
//E.g., for a plane, the light source below the plane cannot cast light on the top, hence, there can also not be any specularity. 
Vec3Df phongSpecularOnly(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, const Vec3Df & cameraPos, unsigned int index)
{
	return Vec3Df(0,1,0);
}

//Blinn-Phong Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
//Be careful!!! The pseudo code does some additional modifications to the formula seen in the course
//Follow the course version and calculate ONLY Ks * pow(dot(N,H), shininess). The definition of H is given on the page above and in the course.
//The same test as before should be used
Vec3Df blinnPhongSpecularOnly(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, const Vec3Df & cameraPos, unsigned int index)
{
	return Vec3Df(0,0,1);
}


///////////////
//Toon Shading
///////////////
//use the variable ToonDiscretize.
//Normal diffuse shading has values between 0 and Kd (the diffuse coefficient).
//In toon shading, these values are discretized.
//This means, there are N (ToonDiscretize) uniform intervals between 0 and Kd - in this example, we assume a single color channel, you should use the values from the vector Kd 
//Let c(0)=0, c(1) ...c(N), c(N+1)=Kd be the boundary values of these intervals.
//For a value v in [c(i), c(i+1)), the function should return (c(i)+c(i+1))/2.
//For v=Kd, return (c(N)+c(N+1))/2, else 0.
Vec3Df toonShadingNoSpecular(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, unsigned int index)
{
	return Vec3Df(0.5,0.5,0);
}

//Toon shading specularity
//The goal is to add white highlights.
//If a channel of Blinn-Phong Specularity has a value bigger or equal to ToonSpecularThreshold, then set it to 1, else to 0.
Vec3Df toonShadingOnlySpecular(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, const Vec3Df & cameraPos, unsigned int index)
{
	return Vec3Df(0,0,1);
}


///////////////
///INTERACTION
///////////////
Vec3Df userInteractionSphere(const Vec3Df & selectedPos, const Vec3Df & camPos)
{
	//RETURN the new light position, defined as follows.
	//selectedPos is a location on the mesh. Use this location to place the light source to cover the location as seen from camPos.
	//Further, the light should be at a distance of 1.5 from the origin of the scene - in other words, located on a sphere of radius 1.5 around the origin.
	return Vec3Df(1,1,1);
}

Vec3Df userInteractionShadow(const Vec3Df & selectedPos, const Vec3Df & selectedNormal, const Vec3Df & lightPos)
{
	//RETURN the new light position such that the light towards the selectedPos is orthgonal to the normal at that location
	//--- in this way, the shading boundary will be exactly at this location.
	//there are several ways to do this, choose one you deem appropriate given the current light position
	//no panic, I will not judge what solution you chose, as long as the above condition is met.
	return Vec3Df(1,0,1);
}

Vec3Df userInteractionSpecular(const Vec3Df & selectedPos, const Vec3Df & selectedNormal, const Vec3Df & lightPos, const Vec3Df & cameraPos)
{
	//RETURN the new light position such that a specularity (highlight) will be located at selectedPos, when viewed from cameraPos and lit from ligthPos.
	//please ensure also that the light is at a distance of 1 from selectedpos! If the camera is on the wrong side of the surface (normal pointing the other way),
	//then just return the original light position.
	//There is only ONE way of doing this!
	return Vec3Df(0,1,1);
}