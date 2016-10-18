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
using namespace std;

vector<vector<unsigned int>> vertexNeighbours;
vector<vector<unsigned int>> edges;
vector<vector<unsigned int>> triangleNeighbours;
vector<vector<pair<unsigned int, bool>>> boundaryEdgeNetwork;
vector<Vec3Df> averages;
vector<Vec3Df> distanceMeshAverage;

 // A boolean that we set to true when we want to draw the baryCenters
bool drawBaryCenters = false;
bool drawNormals = false;
bool drawColoredTriangles = false;
bool drawBoundaries = false;

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

/** Function that gets called on keypress 1 */
void myFunction1(Mesh m) {
	if (drawBaryCenters)
		drawBaryCenters = false;
	else
		drawBaryCenters = true;
	glutPostRedisplay();
}

/** Function that gets called on keypress 2 */
void myFunction2(Mesh m) {
	if (drawNormals)
		drawNormals = false;
	else
		drawNormals = true;
	glutPostRedisplay;
}

Vec3Df calculateSurfaceNormal(Mesh& m, unsigned int triangleIndex) {
	Vec3Df VecU = m.vertices.at(m.triangles.at(triangleIndex).v[0]).p - m.vertices.at(m.triangles.at(triangleIndex).v[1]).p;
	Vec3Df VecV = m.vertices.at(m.triangles.at(triangleIndex).v[1]).p - m.vertices.at(m.triangles.at(triangleIndex).v[2]).p;
	return Vec3Df::crossProduct(VecU, VecV);
}

/** Function that gets called on keypress 3 */
void myFunction3(Mesh m) {
	if (drawColoredTriangles)
		drawColoredTriangles = false;
	else
		drawColoredTriangles = true;
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
		//current triangle = current triangle
		//loop over vertex combinations in current triangle
			//Check matrix. If this combination has value -1: intialise an edge
			//Store the vertices in the edge + the current triangle.
			//Change the value in matrix to the edge index. (change both!)
			
			//Check matrix. If this combination has a value above -1:
			//add the triangle to the existing edge using the edge index in the matrix
	//dont forget to print a 'done'.
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

		for (int edgeIndex = 0; edgeIndex < boundaryEdges.size(); edgeIndex++) {		//Loop over all boundary edges
			if (!boundaryEdges.at(edgeIndex).second) {									//Check if edge is unmarked
				int Enummer = 0;
				unsigned int currentEdge = edgeIndex;									//Make this the current edge
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

void iteratedAveraging(Mesh m, double step) {
	generateVertexNeighbours(m);
	averages.resize(m.vertices.size());
		distanceMeshAverage.resize(m.vertices.size());
	for (int i = 0; i < m.vertices.size(); i++) {
		unsigned int nNeighbours = vertexNeighbours.at(i).size();
		Vec3Df sum = { 0, 0, 0 };
		for (int j = 0; j < nNeighbours; j++) {
			sum += m.vertices.at(vertexNeighbours.at(i).at(j)).p;
		}
		float x = sum[0] / nNeighbours;
		float y = sum[1] / nNeighbours;
		float z = sum[2] / nNeighbours;
		averages.at(i) = {x , y,z };
		distanceMeshAverage.at(i).distance(m.vertices.at(i).p, averages.at(i));
	}

	for (int i = 0; i < m.vertices.size(); i++) {
		m.vertices.at(i).p = m.vertices.at(i).p + step * distanceMeshAverage.at(i);
	}


	//Generate neighbour network 
	//iterations user defined
	//step size
	//for loop
		//Vector<Vec3Df> averages (as big as m.vertices)
		//Vec3Df = average of all neighbours
		//Vector<vec3Df> vector tussen de mesh vertices en de averages
		//
		//
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

void myFunction6(Mesh m) {
	boundaryVertices(m);
	if (drawBoundaries)
		drawBoundaries = false;
	else
		drawBoundaries = true;
}

void myFunctionS(Mesh m) {
	double step = 1;
	iteratedAveraging(m, step);
	glutPostRedisplay();
}

/** Gets called once the mesh has to be drawn.
 Currently calls the meshs draw function, which just
 draws the faces in plain white.
 With your OpenGL knowledge you could write similar
 draw functions that highlight certain vertices,
 edges or faces for better visualization of the
 results of your functions above. */
void draw(Mesh m) {
	 if(!drawColoredTriangles)
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

	if (drawColoredTriangles) {
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

	if (drawBoundaries) {
		glBegin(GL_LINES);
		glColor3f(1.0, 1.0, 0.0);
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


}