// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {
    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    // Create the empty vertex-edge adjacency sparcematrix with size of |E| x |V|
    SparseMatrix<size_t> A0(mesh->nEdges(), mesh->nVertices());
    // Grant access to geometry cache
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    // Assign value to matrix w.r.t the edges info
    for (Edge e : mesh->edges()) {
        // Get the two vertices of the edge
        Vertex v1 = e.firstVertex();
        Vertex v2 = e.secondVertex();

        // Assign the value to the matrix
        A0.insert(geometry->edgeIndices[e], geometry->vertexIndices[v1]) = 1;
        A0.insert(geometry->edgeIndices[e], geometry->vertexIndices[v2]) = 1;
    }

    return A0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    // TODO
    // Create the empty face-edge adjacency sparcematrix with size of |F| x |E|
    SparseMatrix<size_t> A1(mesh->nFaces(), mesh->nEdges());

    // Grant access to geometry cache
    geometry->requireFaceIndices();
    geometry->requireEdgeIndices();

    // Assign value to matrix w.r.t the faces info
    for (Face f : mesh->faces()) {
        // Assign the value to the matrix
        for (Halfedge he : f.adjacentHalfedges()) {
            A1.insert(geometry->faceIndices[f], geometry->edgeIndices[he.edge()]) = 1;
        }
    }

    return A1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    // TODO
    // Construct a vector of length |V|
    Vector<size_t> V = Vector<size_t>::Zero(mesh->nVertices());
    
    // Iterate through the vertices set in MeshSubset
    for (size_t idx : subset.vertices) {
        // Assign the value to the vector
        V(idx) = 1;
    }
    
    return V;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> E = Vector<size_t>::Zero(mesh->nEdges());
    
    for (size_t idx : subset.edges) {
        // Assign the value to the vector
        E(idx) = 1;
    }
    
    return E;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> F = Vector<size_t>::Zero(mesh->nFaces());
    
    for (size_t idx : subset.faces) {
        // Assign the value to the vector
        F(idx) = 1;
    }
    
    return F;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> starVertices(subset.vertices);
    std::set<size_t> starEdges(subset.edges);
    std::set<size_t> starFaces(subset.faces);

    // Add new connected edges in A0 with selected vertices;
    for (size_t idx : starVertices) {
        for (size_t i = 0; i < A0.rows(); i++) {
            if (A0.coeff(i, idx) == 1) {
                starEdges.insert(i);
            }
        }
    }

    // Add connected faces in A1 with selected edges;
    for (size_t idx : starEdges) {
        for (size_t i = 0; i < A1.rows(); i++) {
            if (A1.coeff(i, idx) == 1) {
                starFaces.insert(i);
            }
        }
    }

    MeshSubset starSubset(starVertices, starEdges, starFaces);

    return starSubset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    std::set<size_t> clVertices(subset.vertices);
    std::set<size_t> clEdges(subset.edges);
    std::set<size_t> clFaces(subset.faces);

    for (size_t idx : clFaces) {
        for (size_t i = 0; i < A1.cols(); i++) {
            if (A1.coeff(idx, i) == 1) {
                clEdges.insert(i);
            }
        }
    }

    for (size_t idx : clEdges) {
        for (size_t i = 0; i < A0.cols(); i++) {
            if (A0.coeff(idx, i) == 1) {
                clVertices.insert(i);
            }
        }
    }

    MeshSubset clSubset(clVertices, clEdges, clFaces);
    
    return clSubset;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    // return Cl(St(m)) - St(Cl(m))
    MeshSubset m1 = closure(star(subset.deepCopy()));
    MeshSubset m2 = star(closure(subset.deepCopy()));

    m1.deleteSubset(m2);
    return m1; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return subset.equals(closure(subset.deepCopy())); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (!isComplex(subset)) {
        return -1;
    }

    if (subset.vertices.size() == 0 && subset.edges.size() == 0 && subset.faces.size() == 0) {
        return -1;
    }

    if (subset.vertices.size() > 0 && subset.edges.size() == 0 && subset.faces.size() == 0) {
        return 0;
    }

    // Detect a 2-complex (triangular mash)
    // Each edge should correspond to one face in the subset
    if (subset.faces.size() > 0) {
        for (size_t index : subset.edges) {
            // std::cout << "Checking edge: " << index << std::endl;
            bool found = false;
            for (size_t i : subset.faces) {
                if (A1.coeff(i, index) == 1) {
                    found = true;
                    // std::cout << "Found one: " << i << "," << index << std::endl;
                    break;
                }
            }
            if (!found) {
                // std::cout << "Isolated edge detected" << std::endl;
                return -1;
            }
        }

        // Each vertex should correspond to one edge in the subset
        for (size_t index : subset.vertices) {
            bool found = false;
            for (size_t i : subset.edges) {
                if (A0.coeff(i, index) == 1) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return -1;
            }
        }

        return 2;
    }

    // Detect a 1-complex (connected graph)
    // Each vertex should correspond to one edge in the subset
    for (size_t index : subset.vertices) {
        bool found = false;
        for (size_t i : subset.edges) {
            if (A0.coeff(i, index) == 1) {
                found = true;
                break;
            }
        }
        if (!found) {
            return -1;
        }
    }

    return 1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    int k = isPureComplex(subset);
    if (k == -1) {
        return MeshSubset();
    }

    if (k == 0) {
        return MeshSubset();
    }

    // For 1-complex, delete all vertices that belong to two edges
    if (k == 1) {
        std::set<size_t> bdVertices(subset.vertices);
        std::set<size_t> bdEdges(subset.edges);
        std::set<size_t> bdFaces(subset.faces);

        // Delete all vertices that belong to two edges
        for (size_t idx : subset.vertices) {
            size_t count = 0;
            for (size_t i : subset.edges) {
                if (A0.coeff(i, idx) == 1) {
                    count++;
                }
            }
            if (count == 2) {
                bdVertices.erase(idx);
            }
        }

        bdEdges.clear();
        bdFaces.clear();

        MeshSubset bdSubset(bdVertices, bdEdges, bdFaces);
        return closure(bdSubset);
    }

    // For 2-complex, delete all edges that belong to two faces(together with the vertices)
    // After that, delete all faces
    if (k == 2) {
        std::set<size_t> bdVertices(subset.vertices);
        std::set<size_t> bdEdges(subset.edges);
        std::set<size_t> bdFaces(subset.faces);

        // Delete all edges that belong to two faces
        for (size_t idx : subset.edges) {
            size_t count = 0;
            for (size_t i : subset.faces) {
                if (A1.coeff(i, idx) == 1) {
                    count++;
                }
            }
            if (count == 2) {
                // Remove together with vertices
                for (size_t i = 0; i < A0.cols(); i++) {
                    if (A0.coeff(idx, i) == 1) {
                        bdVertices.erase(i);
                    }
                }
                bdEdges.erase(idx);
            }
        }

        // Delete all faces
        bdFaces.clear();

        MeshSubset bdSubset(bdVertices, bdEdges, bdFaces);
        return closure(bdSubset);
    }

    return MeshSubset(); // placeholder
}