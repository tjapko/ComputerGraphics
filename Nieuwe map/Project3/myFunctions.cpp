/** This file contains five functions that get
 called when the user presses 1-5 and a draw( Mesh m )
 function that gets called everytime the mesh needs to
 be redrawn. 

 All functions get passed a Mesh object, which contains
 the list of vertices and triangles amongst other things.

 In particular, m.vertices is a vector that contains
 a Vertex object for each vertex of the mesh. The vertex
 object simply holds the 3d world coordinates.
 The order of a vertex in this vector determines its
 index.
 m.triangles is a vector that contains one Triangle
 object for each triangle of the mesh, which then
 simply holds three integers, pointing to the indices
 of the vertices that this triangle is made of.

 The draw(Mesh m) function below gets called
 everytime the screen has to be redrawn (e.g.
 when the mouse has been used). For now it calls
 a very simple draw function implemented in the
 Mesh class, but you can change it in order to
 visualize the results of the exercises.
 A redraw can be force by calling
 glutPostRedisplay(), which has to be done
 once your visualization changes.

 These two lists and the draw() function
 should be all you need to solve the
 exercises.

 Feel free to implement your own auxilary functions,
 you are not restricted to the functions shown below.*/

#include "Vec3D.h"
#include "mesh.h"
#include <GL/glut.h>
#include "myFunctions.h"
#include <cmath>
#include <stdlib.h>
using namespace std;

vector<vector<unsigned int>> vertexNeighbours;
vector<vector<unsigned int>> edges;
vector<vector<unsigned int>> triangleNeighbours;
vector<vector<pair<unsigned int, bool>>> boundaryEdgeNetwork;
vector<Vec3Df> averages;
vector<Vec3Df> distanceMeshAverage;
vector<float> smoothingCCVertex;		//smoothing Color Control Vertex
vector<float> smoothingCCTriangle;		//smoothing Color Control Triangle

 // A boolean that we set to true when we want to draw the baryCenters
bool drawBaryCenters = false;
bool drawNormals = false;
bool drawColoredTriangles = false;
bool drawBoundaries = false;
bool smoothBool = false;
bool drawSmoothMap = false;
Mesh smooth;

// NOTE THAT WE PASS THE MESH AS Mesh& m, and not just Mesh m, like I did in the lecture.
// That caused the mesh to be copied everytime this function gets executed and is what
// caused the program to get stuck.
Vec3Df getBaryCenter(Mesh& m, unsigned int triangleIndex) {
	// The triangle with index triangleIndex can be adressed via m.triangles.at(triangleIndex)
	// The index of the first vertex in that triangle is m.triangles.at(triangleIndex).v[0]
	// The position of a vertex with index vertexIndex can be adressed via m.vertices.at(vertexIndex).p
	// This should explain the way we sum up the three vertex coordinates below.
	Vec3Df baryCenter = 1. / 3. * m.vertices.at(m.triangles.at(triangleIndex).v[0]).p
		+ 1. / 3. * m.vertices.at(m.triangles.at(triangleIndex).v[1]).p
		+ 1. / 3. * m.vertices.at(m.triangles.at(triangleIndex).v[2]).p;
	return baryCenter;
}

Vec3Df calculateSurfaceNormal(Mesh& m, unsigned int triangleIndex) {
	Vec3Df VecU = m.vertices.at(m.triangles.at(triangleIndex).v[0]).p - m.vertices.at(m.triangles.at(triangleIndex).v[1]).p;
	Vec3Df VecV = m.vertices.at(m.triangles.at(triangleIndex).v[1]).p - m.vertices.at(m.triangles.at(triangleIndex).v[2]).p;
	return Vec3Df::crossProduct(VecU, VecV);
}

void generateVertexNeighbours(Mesh m) {
	// Change the size of the neighbour array to the size of the number of vertices in the mesh
	vertexNeighbours.resize(m.vertices.size());
	// Loop over all triangles...
	for (unsigned int i = 0; i<m.triangles.size(); i++) {
		Triangle& currentTriangle = m.triangles.at(i);
		// For each vertex within that triangle ...
		for (unsigned int j = 0; j<3; j++) {
			// Store its two neighbours within that triangle:
			// get the index of the current vertex
			int currentVertexIndex = currentTriangle.v[j];
			// and the list of neighbours within the neighbour array
			vector<unsigned int>& currentNeighbours = vertexNeighbours.at(currentVertexIndex);
			for (unsigned int k = 1; k<3; k++) {
				// For each of the two neighbours of the vertex, check if they are already in the list of his neighbours
				if (find(currentNeighbours.begin(), currentNeighbours.end(), currentTriangle.v[(j + k) % 3]) == currentNeighbours.end()) {
					// and if not, add them to this list.
					currentNeighbours.push_back(currentTriangle.v[(j + k) % 3]);
				}
			}
		}
	}
}

void generateTriangleNeighbours(Mesh m) {
	//triangleneigbours initiliase and resize
	triangleNeighbours.resize(m.triangles.size());
	//loop over all triangles
	for (int i = 0; i < m.triangles.size(); i++) {
		Triangle currentTriangle = m.triangles.at(i);
		//loop three times (becasue a triangle can only have three neighbours
		for (int j = 0; j < 3; j++) {
			int vertexA = currentTriangle.v[j];
			int vertexB = currentTriangle.v[(j + 1) % 3];
			for (int k = 0; k < m.triangles.size(); k++) {
				if (i != k) {
					if (m.triangles.at(k).v[0] == vertexA && m.triangles.at(k).v[1] == vertexB) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
					else if (m.triangles.at(k).v[0] == vertexB && m.triangles.at(k).v[1] == vertexA) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
					else if (m.triangles.at(k).v[1] == vertexA && m.triangles.at(k).v[2] == vertexB) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
					else if (m.triangles.at(k).v[1] == vertexB && m.triangles.at(k).v[2] == vertexA) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
					else if (m.triangles.at(k).v[0] == vertexA && m.triangles.at(k).v[2] == vertexB) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
					else if (m.triangles.at(k).v[0] == vertexB && m.triangles.at(k).v[2] == vertexA) {
						triangleNeighbours.at(i).push_back(k);
						k = m.triangles.size();
					}
				}
			}
		}
	}
			// vertex a = currenttriangle[0]
			// same for vertex b
	cout << "Yeeey we did it motherfucker";
			//loop over all triangles
				//if vertex a and b exist in triangle, add it to the neighbours of current triangle
				//if 3 neighbours are found, break out of the for loop;

			
}

void generateEdges(Mesh m) {
	unsigned int someVertex = 0;			//Start with first vertex to make an edge out of.
	generateVertexNeighbours(m);				//Generate neighbour network
	cout << vertexNeighbours.at(5)[2];
	edges.resize(m.vertices.size()*4);
	for (int i = 0; i < edges.size(); i++){
		for (int j = 0; j < 4; j++) {
			edges.at(i).push_back(0);
		}
	}
	int EdgeIndex = 0;

	for (int i = 0; i < m.vertices.size(); i++) {					//Iterate over each vertex
		for (int j = 0; j < vertexNeighbours.at(i).size(); j++) {	//Iterate over each neighbour of this vertex
			unsigned int CurrentNeighbour = vertexNeighbours.at(i)[j];

			//Check if there is an edge in which this combination exists.
			//If there is no edge in which this combination exists, then fill this edge with these values.

			//if it exists, make a boolean true.

			//Make another if statement with this boolean as input
			bool doesItExist = false;
			for (int k = 0; k < EdgeIndex; k++) {//Iterate over all existing edges
				if ((edges.at(k)[0] == i && edges.at(k)[1] == CurrentNeighbour)			//Check if this edge doesn't already exist
					|| (edges.at(k)[1] == i && edges.at(k)[0] == CurrentNeighbour)) {
					doesItExist = true;
					break;
				}
			}
			if(!doesItExist){
				edges.at(EdgeIndex)[0] = i;														//Add current vertex
				edges.at(EdgeIndex)[1] = CurrentNeighbour;								//Add its current neighbour
			//	cout << "Edge " << EdgeIndex << " created with vertices " << i << " and " << CurrentNeighbour << "\n";
				EdgeIndex++;
			}
		}
	}	

	for (int i = 0; i < edges.size(); i++) {									//Iterate over all existing edges
		unsigned int VertexA = edges.at(i)[0];									//Set VertexA to first vertex of edge
		unsigned int VertexB = edges.at(i)[1];									//Set VertexB to second vertex of edge
		for (int j = 0; j < m.triangles.size(); j++) {							//Iterate over all triangles.
			Triangle CurrentTriangle = m.triangles[j];							//Set current triangle.
			if (CurrentTriangle.v[0] == VertexA || CurrentTriangle.v[1] == VertexA || CurrentTriangle.v[2] == VertexA) {		//Check if VertexA is in current triangle
				if (CurrentTriangle.v[0] == VertexB || CurrentTriangle.v[1] == VertexB || CurrentTriangle.v[2] == VertexB) {	//Check if VertexB is in current triangle
					edges.at(i).push_back(j);									//Add the index of the triangle to the edge vector
				}
			}
			if (edges.at(i)[3] != NULL) {										//Quit for-loop if edge vector already has 2 triangles
				j = m.triangles.size();
			}

		}
	}
	cout << "Created edge network";
}

void generateEdgesEfficient(Mesh m) {
	//Generate matrix (vertex * vertex) and fill with -1.
	vector<vector<int>> matrix;
	matrix.resize(m.vertices.size());
	for (int i = 0; i < m.vertices.size(); i++) {
		matrix.at(i).resize(m.vertices.size());
	}
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.at(i).size(); j++) {
			matrix.at(i).at(j) = -1;
		}
	}
	
	edges.resize(3 * m.vertices.size());
	for (int i = 0; i < edges.size(); i++) {
		edges.at(i).resize(4);
		for (int j = 0; j < 4; j++) {
			edges.at(i).at(j) = -1;
		}
	}
	unsigned int edgeIndex = 0;

	//loop over all triangles
	for (int i = 0; i < m.triangles.size(); i++) {
		Triangle currentTriangle = m.triangles.at(i);
		for (int j = 0; j < 3; j++) {
			unsigned int vertexA = currentTriangle.v[j];
			unsigned int vertexB = currentTriangle.v[(j + 1) % 3];

			if (matrix.at(vertexA).at(vertexB) == -1) {
				edges[edgeIndex][0] = vertexA;
				edges[edgeIndex][1] = vertexB;
				edges[edgeIndex][2] = i;
				matrix.at(vertexA).at(vertexB) = edgeIndex;
				matrix.at(vertexB).at(vertexA) = edgeIndex;
				edgeIndex++;
			}

			else {
				edges[matrix.at(vertexA).at(vertexB)][3] = i;
			}
		}
	}

	bool notDone = true;
	while (notDone) {
		if (edges.back().at(0) == 4294967295)
			edges.pop_back();
		else
			notDone = false;
	}
		//current triangle = current triangle
		//loop over vertex combinations in current triangle
			//Check matrix. If this combination has value -1: intialise an edge
			//Store the vertices in the edge + the current triangle.
			//Change the value in matrix to the edge index. (change both!)
			
			//Check matrix. If this combination has a value above -1:
			//add the triangle to the existing edge using the edge index in the matrix
	//dont forget to print a 'done'.
}

void boundaryVerticesOUD(Mesh& m) {
	generateEdgesEfficient(m);					//First generate edges to use in this method

	bool boundaryVertices = false;				
	int amountOfBoundaryEdges = 0;				//Initialise the amount of boundary edges to 0;
	for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++) {			//Loop over all edges
		if (edges.at(edgeIndex)[3] == -1) {										//If this edge has only 1 triangle
			boundaryVertices = true;											//Set boundaryVertices to true.
			amountOfBoundaryEdges++;											//Increment amount of boundary edges
		}
	}

	if (boundaryVertices) {														//Only performs this if there are any boundary vertices.
		boundaryEdgeNetwork.resize(amountOfBoundaryEdges);						//resize the amount of boundary components to the amount of boundary edges. (is too much but ok)
		for (int i = 0; i < boundaryEdgeNetwork.size(); i++) {					//Loop over this network
			boundaryEdgeNetwork.at(i).resize(amountOfBoundaryEdges);			//Resize all vector lists to the amount of boundary edges.
		}

		vector<pair<unsigned int, bool>> boundaryEdges;							//Initialise a temporary list for storing all boundary edges and its markings
		boundaryEdges.resize(amountOfBoundaryEdges);							//resize it
		unsigned int boundaryEdgesIndex = 0;									//Set beginning index
		for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++) {		//Loop over all edges to pick out boundary edges
			if (edges.at(edgeIndex).at(3) == -1) {								// and store them in the temporary list
				boundaryEdges.at(boundaryEdgesIndex).first = edgeIndex;
				boundaryEdges.at(boundaryEdgesIndex).second = false;			// + make their marking false
				boundaryEdgesIndex++;
			}
		}

		int listIndex = 0;				//Index of Boundary components				//Index of Boundary edges per boundary component

		for (int edgeIndex = 0; edgeIndex < boundaryEdges.size(); edgeIndex++) {		//Loop over all boundary edges
			if (!boundaryEdges.at(edgeIndex).second) {									//Check if edge is unmarked
				int Enummer = 0;
				unsigned int currentEdge = boundaryEdges.at(edgeIndex).first;									//Make this the current edge
				boundaryEdgeNetwork[listIndex].at(Enummer).first = currentEdge;			//Add this edge to the final vector
				boundaryEdgeNetwork[listIndex].at(Enummer).second = true;				//Mark 
				boundaryEdges.at(edgeIndex).second = true;
				Enummer++;
				int vertexA = edges.at(boundaryEdges.at(edgeIndex).first).at(0);
				int vertexNext = edges.at(boundaryEdges.at(edgeIndex).first).at(1);
				bool notDone = true;
				while (notDone) {
					for (int i = 0; i < boundaryEdges.size(); i++) {
						if (edges.at(boundaryEdges.at(i).first).at(0) == vertexNext &&
							!boundaryEdgeNetwork[listIndex].at(Enummer - 1).second) {

							boundaryEdgeNetwork[listIndex].at(Enummer).first = boundaryEdges.at(i).first;
							boundaryEdgeNetwork[listIndex].at(Enummer).second = true;
							boundaryEdges.at(Enummer).second = true;
							Enummer++;
							vertexNext = edges.at(boundaryEdges.at(i).first).at(1);
							if (vertexNext == vertexA)
								notDone = false;
						}
						else if (edges.at(boundaryEdges.at(i).first).at(1) == vertexNext &&
							!boundaryEdgeNetwork[listIndex].at(Enummer - 1).second) {

							boundaryEdgeNetwork[listIndex].at(Enummer).first = boundaryEdges.at(i).first;
							boundaryEdgeNetwork[listIndex].at(Enummer).second = true;
							boundaryEdges.at(Enummer).second = true;
							Enummer++;
							vertexNext = edges.at(boundaryEdges.at(i).first).at(0);
							if (vertexNext == vertexA)
								notDone = false;
						}
					}
				}
				listIndex++;
			}
		}
		cout << "The number of boundary components is " << listIndex;

		//Loop over all edges
			//Store boundary edges in seperate list and add a marked int variable

		//Loop over all boundary edges
			// check if edge is unmarked
				//make edge 'current edge' and mark it.
				//store this edge in list
				//Store first vertex of this edge in memory and use second vertex to find next boundary edge
				//Loop over all boundary edges
					//Find the boundary edge next to current edge, until first vertex of first boundary edge is found. (Should not find itself!)
					//mark the next edge
		
		//Output the size of the list<lists>

		

		//draw this shit


	}
	
}

void boundaryVertices(Mesh& m) {
	generateEdgesEfficient(m);					//First generate edges to use in this method

	bool boundaryVertices = false;
	int amountOfBoundaryEdges = 0;				//Initialise the amount of boundary edges to 0;
	for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++) {			//Loop over all edges
		if (edges.at(edgeIndex)[3] == -1) {										//If this edge has only 1 triangle
			boundaryVertices = true;											//Set boundaryVertices to true.
			amountOfBoundaryEdges++;											//Increment amount of boundary edges
		}
	}

	if (boundaryVertices) {														//Only performs this if there are any boundary vertices.
		boundaryEdgeNetwork.resize(amountOfBoundaryEdges);						//resize the amount of boundary components to the amount of boundary edges. (is too much but ok)
		for (int i = 0; i < boundaryEdgeNetwork.size(); i++) {					//Loop over this network
			boundaryEdgeNetwork.at(i).resize(amountOfBoundaryEdges);			//Resize all vector lists to the amount of boundary edges.
		}

		vector<pair<unsigned int, bool>> boundaryEdges;							//Initialise a temporary list for storing all boundary edges and its markings
		boundaryEdges.resize(amountOfBoundaryEdges);							//resize it
		unsigned int boundaryEdgesIndex = 0;									//Set beginning index
		for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++) {		//Loop over all edges to pick out boundary edges
			if (edges.at(edgeIndex).at(3) == -1) {								// and store them in the temporary list
				boundaryEdges.at(boundaryEdgesIndex).first = edgeIndex;
				boundaryEdges.at(boundaryEdgesIndex).second = false;			// + make their marking false
				boundaryEdgesIndex++;
			}
		}

		int listIndex = 0;				//Index of Boundary components				//Index of Boundary edges per boundary component

		for (int edgeIndex = 0; edgeIndex < boundaryEdges.size(); edgeIndex++) {		//Loop over all boundary edges through temporary BoundaryEdges list
			if (!boundaryEdges.at(edgeIndex).second) {									//Check if edge is unmarked
				int Enummer = 0;
				unsigned int currentEdge = boundaryEdges.at(edgeIndex).first;			//Make this the current edge
				boundaryEdgeNetwork[listIndex].at(Enummer).first = currentEdge;			//Add this edge to the final vector list
				boundaryEdgeNetwork[listIndex].at(Enummer).second = true;				//Mark this edge
				boundaryEdges.at(edgeIndex).second = true;								//Mark this edge in temporary list
				Enummer++;																//Increment Enummer to find next edge
				int vertexA = edges.at(currentEdge).at(0);								//Store first vertex of current edge
				int vertexNext = edges.at(currentEdge).at(1);							//Use second vertex to find the next edge
				bool notDone = true;
				while (notDone) {
					for (int i = 0; i < boundaryEdges.size(); i++) {					//Loop over all boundary edges
						if (edges.at(boundaryEdges.at(i).first).at(0) == vertexNext &&		//Checks if this edge consists of the second vertex of previous edge
							boundaryEdgeNetwork[listIndex].at(Enummer - 1).first != boundaryEdges.at(i).first) {	// checks if we do not find the same edge

							boundaryEdgeNetwork[listIndex].at(Enummer).first = boundaryEdges.at(i).first;		//Add this edge index to the boundary edges vector
							boundaryEdgeNetwork[listIndex].at(Enummer).second = true;							//Mark this edge
							boundaryEdges.at(Enummer).second = true;											//Mark this edge
							Enummer++;																			//Increment Enummer to find next edge
							vertexNext = edges.at(boundaryEdges.at(i).first).at(1);								//Set vertexNext to find the next edge
							if (vertexNext == vertexA)															//If the first vertex is found, break out of while loop
								notDone = false;
						}
						else if (edges.at(boundaryEdges.at(i).first).at(1) == vertexNext &&			//Same as above, but checks for the second vertex of next edge
							boundaryEdgeNetwork[listIndex].at(Enummer - 1).first != boundaryEdges.at(i).first) {

							boundaryEdgeNetwork[listIndex].at(Enummer).first = boundaryEdges.at(i).first;
							boundaryEdgeNetwork[listIndex].at(Enummer).second = true;
							boundaryEdges.at(Enummer).second = true;
							Enummer++;
							vertexNext = edges.at(boundaryEdges.at(i).first).at(0);					//Sets the first vertex as next vertex
							if (vertexNext == vertexA)
								notDone = false;
						}
					}
				}
				listIndex++;												//When the boundary edge component is round, increment listindex
			}
		}
		cout << "The number of boundary components is " << listIndex << "\n";

	}

}

int connectedComponents(Mesh& m) {
	unsigned int res;					//Set result.
	res = 0;

	vector<unsigned int> processMe;						//Create a Process vector
	vector<bool> markedVertices;						//Create another vector that contains all the marked/processed vertices
	markedVertices.resize(m.vertices.size());			//Make it as big as the total number of vertices
	for (int i = 0; i < markedVertices.size(); i++)		//set every vertex to marked = false.
		markedVertices[i] = false;

	for (int i = 0; i < markedVertices.size(); i++) {	//Iterate over all vertices
		if (!markedVertices[i]) {						//If this vertex is not marked, add it to the processMe vector
			processMe.push_back(i); 
			while (processMe.size() > 0) {									//while there are still vertices you need to process
				unsigned int nextVertexToProcess = processMe.back();		//gets the last vertex in the list
				if (markedVertices[nextVertexToProcess] == false) markedVertices[nextVertexToProcess] = true;	//Mark the vertex
				processMe.pop_back();										//removes the last element of this list

																									//Now the neighbours need to be added.
				for (int j = 0; j < vertexNeighbours.at(nextVertexToProcess).size(); j++) {			//Iterate over all neighbours
					unsigned int currentNeighbour = vertexNeighbours.at(nextVertexToProcess)[j];	//Set current Neighbour
					if (markedVertices[currentNeighbour] == false)	//if this neighbour is already marked, discard.
						processMe.push_back(currentNeighbour);		//otherwise, add it to processMe.
				}
			}
			res = res + 1;
		}
	}
	return res;
}

Mesh smoothing(Mesh m, double step) {
	Mesh newMesh = m;
	generateVertexNeighbours(newMesh);
	averages.resize(newMesh.vertices.size());
	distanceMeshAverage.resize(newMesh.vertices.size());
	for (int i = 0; i < newMesh.vertices.size(); i++) {
		unsigned int nNeighbours = vertexNeighbours.at(i).size();
		Vec3Df sum = { 0, 0, 0 };
		for (int j = 0; j < nNeighbours; j++) {
			sum += newMesh.vertices.at(vertexNeighbours.at(i).at(j)).p;
		}
		float x = sum[0] / nNeighbours;
		float y = sum[1] / nNeighbours;
		float z = sum[2] / nNeighbours;
		averages.at(i) = {x , y,z };
		distanceMeshAverage.at(i) = averages.at(i) - newMesh.vertices.at(i).p;
	}

	if (step > 0) {
		for (int i = 0; i < newMesh.vertices.size(); i++) {
			newMesh.vertices.at(i).p = newMesh.vertices.at(i).p + step * distanceMeshAverage.at(i);
		}
	}
	

	//Generate neighbour network 
	//iterations user defined
	//step size
	//for loop
		//Vector<Vec3Df> averages (as big as m.vertices)
		//Vec3Df = average of all neighbours
		//Vector<vec3Df> vector tussen de mesh vertices en de averages
	return newMesh;
}

float getMin(vector<float> A) {
	float B = A.at(0);
	for (int i = 0; i < A.size(); i++) {
		if (A.at(i) < B)
			B = A.at(i);
	}
	return B;
}

float getMax(vector<float> A) {
	float B = A.at(0);
	for (int i = 0; i < A.size(); i++) {
		if (A.at(i) > B)
			B = A.at(i);
	}
	return B;
}

void generateSmoothMap(Mesh m) {
	smoothingCCVertex.resize(m.vertices.size());
	for (int i = 0; i < smoothingCCVertex.size(); i++) {
		smoothingCCVertex.at(i) = distanceMeshAverage.at(i).getLength();
	}

	float min = getMin(smoothingCCVertex);
	float max = getMax(smoothingCCVertex);

	for (int i = 0; i < smoothingCCVertex.size(); i++) {
		smoothingCCVertex.at(i) = (smoothingCCVertex.at(i) - min) / (max - min);
	}

	smoothingCCTriangle.resize(m.triangles.size());
	for (int i = 0; i < smoothingCCTriangle.size(); i++) {
		float sum = 0.0;
		for (int j = 0; j < 3; j++) {
			sum += smoothingCCVertex.at(m.triangles.at(i).v[j]);
		}
		smoothingCCTriangle.at(i) = sum / 3.0;
	}
}

float generateVolume(Mesh m) {
	// pick a point p in space ( PICK THE ORIGEN!!!!!!!!)
	// create a vector triangleVolumes with the size of all triangles
	vector<float> triangleVolumes;
	triangleVolumes.resize(m.triangles.size());
	float totalVolume = 0;
	// loop over all triangles. 
	for (int i = 0; i < triangleVolumes.size(); i++) {
		Vec3Df a = m.vertices.at(m.triangles.at(i).v[0]).p;
		Vec3Df b = m.vertices.at(m.triangles.at(i).v[1]).p;
		Vec3Df c = m.vertices.at(m.triangles.at(i).v[2]).p;
		Vec3Df cross = Vec3Df::crossProduct(b, c);
		float dot = Vec3Df::dotProduct(a, cross);
		float abs = fabs(dot);

		Vec3Df normal = calculateSurfaceNormal(m, i);
		float temp = Vec3Df::dotProduct(normal, a);
		if (temp < 0) 
			abs = abs * -1.0;

		totalVolume += abs / 6.0;
	}
	return totalVolume;
	// for each triangle compute the volume in regards to point p if p is origin is V= 1/6*|a*(b x c)|		a dot (b cross c)
	// and assign and + or - dependeing on if the normal of the surface points towards point p or not. 
	// you do this by taking one of the three verticies and the normal and use the dotproduct between them.  if the outcome of that is positive make the volume negative.
	// put the volmune with the sign in the vector
	// volume it the sum of all values in the vector trianlgeVolumes 
}

Mesh inflate(Mesh m, float initialVolume) {
	//calculate current volume
	//smoothing
	//start while loop
	Mesh tempMesh = m;
	float currentVolume = generateVolume(m);
	float step = 1.001;
	int iteration = 1;
	while (currentVolume < initialVolume) {
		for (int i = 0; i < tempMesh.vertices.size(); i++) {
			tempMesh.vertices.at(i).p[0] = tempMesh.vertices.at(i).p[0] * step;
			tempMesh.vertices.at(i).p[1] = tempMesh.vertices.at(i).p[1] * step;
			tempMesh.vertices.at(i).p[2] = tempMesh.vertices.at(i).p[2] * step;
		}
		currentVolume = generateVolume(tempMesh);
		iteration++;
	}
	cout << "Volume after " << iteration << " iterations with step-size " << step << ": " << currentVolume << "\n";
	cout << "The overshoot in volume is " << (currentVolume - initialVolume) << "\n";
	cout << "Inflating done\n";
	return tempMesh;
		//Use step to increase vertex positions
		//Check if volume is back to normal
			//if it is, break while loop
	//output mesh

}

/** Function that gets called on keypress h */
void myFunctionH() {
	cout << "What would you like to do:\n1 - Toggle barycenters\n2 - Toggle normals\n3 - Toggle colored triangles\n4 - Generate VertexNeighbours and show number of connected components\n" <<
		"5 - Generate Edge Network and show number of edges\n6 - Toggle boundary edges and show number of connected boundary components\n" <<
		"7 - Toggle SmoothMap\n8 - Generate and output volume of the mesh\n9 - Smooth and inflate mesh\nw - Toggle triangles\ns - Smooth the mesh\nh - help\n";
}

/** Function that gets called on keypress 1 */
void myFunction1(Mesh m) {
	if (drawBaryCenters) {
		drawBaryCenters = false;
		cout << "drawBaryCenters OFF\n";
	}
	else {
		drawBaryCenters = true;
		cout << "drawBaryCenters ON\n";
	}
	glutPostRedisplay();
}

/** Function that gets called on keypress 2 */
void myFunction2(Mesh m) {
	if (drawNormals) {
		drawNormals = false;
		cout << "drawNormals OFF\n";
	}
	else {
		drawNormals = true;
		cout << "drawNormals ON\n";
	}
	glutPostRedisplay;
}

/** Function that gets called on keypress 3 */
void myFunction3(Mesh m) {
	if (drawColoredTriangles) {
		drawColoredTriangles = false;
		cout << "drawColoredTriangles OFF\n";
	}
	else {
		drawColoredTriangles = true;
		cout << "drawColoredTriangles ON\n";
	}
}

/** Function that gets called on keypress 4 */
void myFunction4(Mesh m) {
	generateVertexNeighbours(m);
	cout << "This mesh consists of " << connectedComponents(m) << " connected component(s)\n";
}

/** Function that gets called on keypress 5 */
void myFunction5(Mesh m) {
	generateEdgesEfficient(m);
	cout << "Edge Network created succesfully. This mesh has " << edges.size() << " edges.\n";
}

/** Function that gets called on keypress 6 */
void myFunction6(Mesh m) {
	boundaryVertices(m);
	if (drawBoundaries) {
		drawBoundaries = false;
		cout << "drawBoundaries OFF\n";
	}
	else {
		drawBoundaries = true;
		cout << "drawBoundaries ON\n";
	}
}

/** Function that gets called on keypress 7 */
void myFunction7(Mesh m) {
	generateSmoothMap(m);
	if (drawSmoothMap) {
		drawSmoothMap = false;
		cout << "SmoothMap OFF\n";
	}
	else {
		drawSmoothMap = true;
		cout << "SmoothMap ON\n";
	}
}

/** Function that gets called on keypress 8 */
void myFunction8(Mesh m) {
	float volume = generateVolume(m);
	cout << "The volume of this mesh is " << volume << "\n";
}

/** Function that gets called on keypress 9 */
Mesh myFunction9(Mesh m) {
	float initialVolume = generateVolume(m);
	cout << "The initial volume is " << initialVolume << "\n";
	smooth = myFunctionS(m);
	smooth = inflate(smooth, initialVolume);
	draw(smooth);
	return smooth;
}

/** Function that gets called on keypress 0 */
void myFunction0(Mesh m) {}

/** Function that gets called on keypress s */
Mesh myFunctionS(Mesh m) {
	cout << "Define a step-size for the averaging iteration (Enter a number)\n";
	double step;
	cin >> step;
	cout << "You entered step-size: " << step << "\n";
	cout << "How many averaging iterations would you like to perform? (Enter a number)\n";
	int X;
	cin >> X;
	cout << "You entered " << X << " iterations\n";
	for (int i = 0; i < X; i++) {
		smooth = smoothing(m, step);
		draw(smooth);
		glutPostRedisplay();
		smoothBool = true;
	}
	cout << "Smoothing complete\n";
	return smooth;
}

/** Gets called once the mesh has to be drawn.
 Currently calls the meshs draw function, which just
 draws the faces in plain white.
 With your OpenGL knowledge you could write similar
 draw functions that highlight certain vertices,
 edges or faces for better visualization of the
 results of your functions above. */
void draw(Mesh m) {
	if (smoothBool) {
		m = smooth;
	}

	if(!drawColoredTriangles && !drawBoundaries && !drawSmoothMap)
		 m.draw();

	if (drawBaryCenters) {
		// Tell OpenGL to draw large dots
		// WARNING: This has to appear before glBegin()!
		// you cannot vary point size in between one glBegin()/glEnd() block!
		glPointSize(6.);
		// Tell OpenGL to draw poitns
		glBegin(GL_POINTS);
		// Tell OpenGL to draw red dots
		glColor3f(1., 0, 0);
		// Draw points at the position of the barycenters:
		for (int i = 0; i<m.triangles.size(); i++) {
			Vec3Df bary = getBaryCenter(m, i);
			glVertex3f(bary[0], bary[1], bary[2]);
		}
		glEnd();
		// Start drawing in white again to prevent the mesh from being red
		glColor3f(1., 1., 1.);
	}

	if (drawNormals) {
		glBegin(GL_LINES);
		glColor3f(1.0, 1.0, 0.0);
		for (int i = 0; i < m.triangles.size(); i++) {
			Vec3Df bary = getBaryCenter(m, i);
			Vec3Df norm = calculateSurfaceNormal(m, i);
			glVertex3f(bary[0], bary[1], bary[2]);
			glVertex3f((bary[0]+15*norm[0]), (bary[1]+15*norm[1]), (bary[2]+15*norm[2]));
		}
		glEnd();
		glColor3f(1., 1., 1.);

	}

	//DrawColoredTriangles oud
/*	if (drawColoredTriangles3) {
		glBegin(GL_TRIANGLES);
		float x = 1.0 / m.triangles.size();
		for (int i = 0; i < m.triangles.size(); i++) {
			glColor3f(0.0 + x*i, 0.5, 0.5);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}
		glEnd();
		glColor3f(1.0, 1.0, 1.0);
	}
	*/

	if (drawColoredTriangles) {
		unsigned int n = m.triangles.size();
		float nPart = n / 6.0;
		glBegin(GL_TRIANGLES);

		//Part 1 of 6
		for (int i = 0; i < nPart; i++) {
			glColor3f(1.0, 0.0, 1.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}

		//Part 2 of 6
		for (int i = nPart; i < 2*nPart; i++) {
			glColor3f(0.0, 0.0, 1.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}

		//Part 3 of 6
		for (int i = 2*nPart; i < 3 * nPart; i++) {
			glColor3f(0.0, 1.0, 1.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}

		//Part 4 of 6
		for (int i = 3 * nPart; i < 4 * nPart; i++) {
			glColor3f(0.0, 1.0, 0.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}

		//Part 5 of 6
		for (int i = 4 * nPart; i < 5 * nPart; i++) {
			glColor3f(1.0, 1.0, 0.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}

		//Part 6 of 6
		for (int i = 5 * nPart; i < n; i++) {
			glColor3f(1.0, 0.0, 0.0);
			float Xtriangle = m.vertices.at(m.triangles.at(i).v[0]).p[0];
			float Ytriangle = m.vertices.at(m.triangles.at(i).v[0]).p[1];
			float Ztriangle = m.vertices.at(m.triangles.at(i).v[0]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[1]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[1]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[1]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);

			Xtriangle = m.vertices.at(m.triangles.at(i).v[2]).p[0];
			Ytriangle = m.vertices.at(m.triangles.at(i).v[2]).p[1];
			Ztriangle = m.vertices.at(m.triangles.at(i).v[2]).p[2];
			glVertex3f(Xtriangle, Ytriangle, Ztriangle);
		}
		glEnd();
		glColor3f(1.0, 1.0, 1.0);
	}

	if (drawBoundaries) {
		m.draw();
		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);
		for (int j = 0; j < boundaryEdgeNetwork.size(); j++) {
			for (int i = 0; i < boundaryEdgeNetwork.at(j).size(); i++) {
				Vec3Df vertexA = m.vertices.at(edges.at(boundaryEdgeNetwork.at(j).at(i).first).at(0)).p;
				Vec3Df vertexB = m.vertices.at(edges.at(boundaryEdgeNetwork.at(j).at(i).first).at(1)).p;
				glVertex3f(vertexA[0], vertexA[1], vertexA[2]);
				glVertex3f(vertexB[0], vertexB[1], vertexB[2]);
			}
		}
		glEnd();
		glColor3f(1., 1., 1.);
	}

	if (drawSmoothMap) {
		glBegin(GL_TRIANGLES);
		for (int triangleIndex = 0; triangleIndex < smooth.triangles.size(); triangleIndex++) {
			glColor3f(0 + (smoothingCCTriangle.at(triangleIndex) * 2.0), 2.0 - (smoothingCCTriangle.at(triangleIndex) * 2.0), 0.0);
			for (int i = 0; i < 3; i++) {
				glVertex3f(smooth.vertices.at(smooth.triangles.at(triangleIndex).v[i]).p[0],
					smooth.vertices.at(smooth.triangles.at(triangleIndex).v[i]).p[1],
					smooth.vertices.at(smooth.triangles.at(triangleIndex).v[i]).p[2]);
			}
		}
		glEnd();
		glColor3f(1., 1., 1.);
	}

}