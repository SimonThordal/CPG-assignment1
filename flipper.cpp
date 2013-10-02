/* GEL demo program.
 
 Creates a very small simple mesh as a Manifold (cf. the HMesh namespace of the 
 GEL library). The mesh is visualized and pressing any key will cause an edge to
 flip (if the edge is legally flippable). While the example seems pointless at first
 glance it provides a good starting point for implementing 2D Delaunay triangulation
 and also includes point loading code. 
 
 Andreas Bærentzen and Kasper Steenstrup 2013
 
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <HMesh/Manifold.h>
#include <HMesh/AttributeVector.h>
#include <HMesh/triangulate.h>

#include <GLGraphics/gel_glut.h>
#include <CGLA/Vec2d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat3x3f.h>
#include <CGLA/Mat4x4f.h>

using namespace std;
using namespace CGLA;
using namespace HMesh;

// The range of the input data.
Vec2d dmin(99e99), dmax(-99e99);

Manifold m; // The triangle mesh data structure.
HalfEdgeIDIterator flipper = m.halfedges_begin(); // The halfedge we try to flip.
HalfEdgeAttributeVector<int> touched;

/*
 * Draw the triangle mesh. This function is called from GLUT (a library
 * which enables OpenGL drawing in windows.)
 */

void display()
{
	// Set up correct OpenGL projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(dmin[0], dmax[0], dmin[1], dmax[1]);
	glMatrixMode(GL_MODELVIEW);

	// Specify that we want to draw triangle outlines
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Black on white.
	glClearColor(1,1,1,0);
	glColor3f(0,0,0);

	// Clear the screen.
	glClear(GL_COLOR_BUFFER_BIT);

  for(FaceID f: m.faces()){
    glBegin(GL_POLYGON);
    for(Walker w = m.walker(f); !w.full_circle(); w = w.next())
      glVertex3dv(m.pos(w.vertex()).get());
    glEnd();
  }

	// Draw flipper.
	glColor3f(1,0,0);
	glBegin(GL_LINES);
  Walker hew = m.walker(*flipper);
  glVertex3dv(m.pos(hew.vertex()).get());
	glVertex3dv(m.pos(hew.opp().vertex()).get());
	glEnd();

	glFinish();
}

/*
 * As the name indicates, the function below creates a manifold (i.e. a 
 * triangle mesh) consisting of just a single triangle.
 */

void create_single_triangle_manifold(const Vec3f& p1, 
																		 const Vec3f& p2, 
																		 const Vec3f& p3, 
																		 Manifold& mani)
{
	// Create vector of vertices
	vector<Vec3f> vertices(3);
	vertices[0] = p1;
	vertices[1] = p2;
	vertices[2] = p3;

	// Create vector of faces. Each element corresponds to a face and tells
	// how many vertices that face contains. In the case of a triangle
	// mesh, each face has three vertices.
	vector<int> faces(1);
	faces[0] = 3;

	// Create the index vector. Each element is an index into the vertex list
	// 
	vector<int> indices(3);
	indices[0]=0;
	indices[1]=1;
	indices[2]=2;

  mani.build(3,           // Number of vertices.
	  				 reinterpret_cast<float*>(&vertices[0]),// Pointer to vertices.
						 1,           // Number of faces.
						 &faces[0],   // Pointer to faces.
						 &indices[0]);// Pointer to indices.


}


void keyfun(unsigned char c, int x, int y)
{
	/*
	 * A little game, try to flip an edge when user presses any key.
	 * Not all edges can be flipped. Boundary edges cannot. Edges also
	 * cannot be flipped if it will render the mesh invalid.
	 */

	if(boundary(m, *flipper)) // If this is a boundary edge just drop the idea.
		cout << "boundary edge" << endl;
	else if(precond_flip_edge(m, *flipper))
  {
    m.flip_edge(*flipper);
		cout << "flipped" << endl;
  }
	else
  	cout << "could not flip" << endl;

	do
		{
			++flipper; // Get the next halfedge
			// If we have passed the last halfedge, go to the first.
			if(flipper==m.halfedges_end())
				{
					flipper = m.halfedges_begin();
					break;
				}
		}
	while(touched[*flipper] == 0); // Only visit halfedges marked '1'


	// Function call below informs glut that display should be called to
	// show the window again.
	glutPostRedisplay();
}


/*
 * left of
 */
bool leftOf(const Vec3d p1, const Vec3d p2, const Vec3d testPoint) {
	Mat3x3f lo( Vec3f(1), Vec3f(p1[0], p2[0], testPoint[0]), Vec3f(p1[1], p2[1], testPoint[1]) );
	// Find determinanten
	float d = determinant(lo);
	return d > 0;
}

/*
* Check if a point, p4, is inside the circumcirle of p1, p2 and p3
*/
bool inCircle(const Vec3d p1, const Vec3d p2, const Vec3d p3, const Vec3d p4) {
	Mat4x4f pointMatrix( Vec4f(1,1,1,1), Vec4f(p1[0], p2[0], p3[0], p4[0]), Vec4f(p1[1], p2[1], p3[1], p4[1]), Vec4f( pow(p1[0], 2) + pow(p1[1], 2), pow(p2[0], 2) + pow(p2[1], 2), pow(p3[0], 2) + pow(p3[1], 2), pow(p4[0], 2) + pow(p4[1], 2)) );
	// Find determinant
	float d = determinant(pointMatrix);
	return d < 0;
}

/*
* Checks if an edge is a boundary and if not, checks it for being locally Delaunay.
* If it is not locally delaunay it is flipped and all its neighbor edges are checked, ad nauseam.
*/
void recursiveDelaunayFlip(Manifold &m, Walker w, bool isAffected) {
	if ( !boundary(m, w.halfedge()) ) {
		// Check if the current halfedge is locally Delaunay using the inCircle function
		Vec3d p1 = m.pos(w.opp().vertex());
		Vec3d p2 = m.pos(w.vertex());
		Vec3d p3 = m.pos(w.next().vertex());
		Vec3d p4 = m.pos(w.opp().next().vertex()); // This seems to return erroneus values every time
		if (isAffected == true) {
			cout << "Affected quadrillateral:" << endl;
			cout << "p1: " << p1 << endl;
			cout << "p2: " << p2 << endl;
			cout << "p3: " << p3 << endl;
			cout << "p4: " << p4 << endl;
		}
		if ( inCircle( p1, p3, p2, p4 ) || inCircle(p1, p2, p4, p3) ) {
			// Since either point was in a triangle circumcircle, flip the edge
			m.flip_edge(w.halfedge());
			cout << "Edge to be flipped: " << p1 << ", " << p2 << ". Other vertices: " << p3 << ", " << p4 << endl;
			// Recursively check all the edges that share a neighbour with the flipped edge
			recursiveDelaunayFlip(m, w.next(), true);
			recursiveDelaunayFlip(m, w.prev(), true);
			recursiveDelaunayFlip(m, w.opp().next(), true);
			recursiveDelaunayFlip(m, w.opp().prev(), true);
		}
	}
}

/*
 * mark halfedges.
 *
 * Remember that a geometric edge corresponds to two halfedges.
 * this function loops over all halfedges and tags precisely
 * one halfedge in every pair of halfedges with '1'. Its opposite edge
 * is tagged with '0'.
 */

void mark_halfedges()
{
	// Give all halfedges a mark of 0
  for(HalfEdgeID h: m.halfedges())
  {
	  if(m.walker(h).opp().halfedge() < h)
      touched[h] = 0;
    else
      touched[h] = 1;
    }
}


int main(int argc, char** argv)
{
	/*
	 * Read and parse a point set.
	 */

	/* Open a data stream for reading.
	 * We first open data.txt. There is also kote1.txt which contains height
	 * values in addition to x,y positions.
	 */

	ifstream data("data.txt");
	vector<Vec2d> pts;
	if(data.good())
		while(!data.eof())
			{
				double x,y;
				data >> x >> y;

				if(data.good())
					{
						Vec2d p(x,y);
						pts.push_back(p);
						dmin = v_min(p,dmin);
						dmax = v_max(p,dmax);
					}
			}
	cout << "Loaded " << pts.size() << " points "  <<  endl;

	Vec2d trans((dmax[0]+dmin[0])/2,(dmax[1]+dmin[1])/2);
	double skal = 2/max(dmax[0]-dmin[0],dmax[1]-dmin[1]);

	/* Træk trans fra alle punkter og gang med 'skal'*/
	for (int i = 0; i < pts.size(); i++) {
		pts[i] -= trans;
		pts[i] *= skal;
	}

	/*
	 * Build a triangle mesh with a single triangle consisting of the
	 * first three vertices.
	 */

	create_single_triangle_manifold(Vec3f(0, 3, 0),
																	Vec3f(4.5, -1.5, 0),
																	Vec3f(-4.5, -1.5, 0),
																	m);
	// Initially just split the triangle by inserting the first point
	VertexID v = m.split_face_by_vertex(*m.faces_begin());
	m.pos(v) = Vec3d(pts[0][0], pts[0][1], 0);
	// Now insert all of the remaining points
	for (int i = 1; i < pts.size(); i++) {		
		Vec3d insertionPoint = Vec3d(pts[i][0], pts[i][1], 0);
		VertexID insertionVertex;

		// Loop over all the faces and find the face that contains the point
		for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); f++) {
			Walker w = m.walker(*f);
			bool isLeftOf = true;
			while (!w.full_circle()) {
				// If the point to be inserted is not to the left of the halfedge, then break the while loop and continue to the next face
				if (!leftOf(m.pos(w.circulate_face_ccw().vertex()), m.pos(w.vertex()), insertionPoint)) {
					isLeftOf = false;
					break;
				}
				w = w.circulate_face_cw();
			}

			// if we found the face the point belongs to then insert it and break the for loop
			if (isLeftOf == true) {
				insertionVertex = m.split_face_by_vertex(*f);
				m.pos(insertionVertex) = insertionPoint;
				break;
			}
		}

		// Now loop over all the edges affected by the inserted point. 
		// Note that we are assuming that the point was inserted, if not then this will crash spectacularly.
		Walker w = m.walker(insertionVertex);
		// Keep track of the next halfedge pointing TO the inserted vertex
		HalfEdgeID next_edge = w.circulate_vertex_ccw().opp().halfedge();
		HalfEdgeAttributeVector<int> touched;
		while (!w.full_circle()) {
			// Iterate over the face of the current halfedge until we reach the next edge pointing TO the inserted vertex
			if(w.halfedge() != next_edge) {
				// Check if the current halfedge is locally Delaunay using the inCircle function
				recursiveDelaunayFlip(m, w, false);
				// Update the walker to be the next halfedge in the current face.
				w = w.circulate_face_ccw();
			} else {
				// If we are the next edge pointing to the inserted vertex then go to opposite halfedge. This means we are now looking at the halfedge pointing AWAY from the inserted vertex.
				w = w.opp();
				// Remember to update the next_edge to be the next halfedge pointing to the inserted vertex.
				next_edge = w.circulate_vertex_ccw().opp().halfedge();
			}
		}

	}
	/*
	 * Initialize GLUT, the system used to show OpenGL windows.
	 */
	glutInit(&argc, argv);
	glutInitWindowSize(512,512);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow("Delaunay");

	glutDisplayFunc(display); // This function is called from glut to draw
	glutKeyboardFunc(keyfun); // Parse keyboard input
	
	// Pass control to glut
	glutMainLoop();
	return 0;
}
